#define USE_FC_LEN_T
#include <algorithm>
#include <cmath>
#include <cfloat>
#include "psis.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#ifndef FCONE
# define FCONE
#endif

/*****************************************
 Pareto smoothed importance sampling (PSIS) leave-one-out predictive densities.

 Follows the algorithm of Vehtari, Simpson, Gelman, Yao and Gabry (2024),
 "Pareto smoothed importance sampling", JMLR 25(72), as implemented in the R
 packages loo (>= 2.10.1; psis(), loo()) and posterior (>= 1.7.0; gpdfit(),
 qgeneralized_pareto()), with relative efficiency r_eff = 1 (independent
 posterior draws). Each observation is processed with O(S) workspace supplied
 by the caller; no n x S matrix is formed.
 *****************************************/

// long double accumulation, as in base R sum()/mean() and matrixStats
#define LDOUBLE long double

// Tail length used to fit the generalized Pareto distribution;
// loo:::n_pareto() with r_eff = 1: ceiling(min(0.2 * S, 3 * sqrt(S)))
int psis_tail_length(int S){
  return (int) ceil(fmin2(0.2 * S, 3.0 * sqrt((double) S)));
}

// Number of grid points in gpdfit(): min_grid_pts + floor(sqrt(N)), min_grid_pts = 30
int psis_gpd_grid_length(int L){
  return 30 + (int) floor(sqrt((double) L));
}

// log(sum(exp(a))), as in matrixStats::logSumExp()
static double psis_logSumExp(double *a, int n){
  int i;
  double a_max = R_NegInf;
  for(i = 0; i < n; i++){
    if(a[i] > a_max) a_max = a[i];
  }
  if(!R_FINITE(a_max)) return a_max;
  LDOUBLE sum = 0.0;
  for(i = 0; i < n; i++){
    sum += exp(a[i] - a_max);
  }
  return a_max + (double) log((double) sum);
}

// mean(x) with the two-pass refinement of base R mean.default()
static double psis_mean(double *x, int n){
  int i;
  LDOUBLE s = 0.0;
  for(i = 0; i < n; i++) s += x[i];
  s /= n;
  if(R_FINITE((double) s)){
    LDOUBLE t = 0.0;
    for(i = 0; i < n; i++) t += (x[i] - s);
    s += t / n;
  }
  return (double) s;
}

/*****************************************
 Fit the generalized Pareto distribution (location 0) to a sample x of size N
 sorted in ascending order, following Zhang and Stephens (2009) with the weakly
 informative prior on k (Vehtari et al. 2024); mirrors posterior::gpdfit(x,
 wip = TRUE, min_grid_pts = 30, sort_x = FALSE). Here k is the negative of k in
 Zhang and Stephens (2009). On exit, x is overwritten (used as workspace).
 theta and l_theta are workspace of length psis_gpd_grid_length(N).
 Returns k = NA, sigma = NA if the first quartile equals the minimum.
 *****************************************/
void psis_gpdfit(double *x, int N, double *theta, double *l_theta,
                 double *k_out, double *sigma_out){

  int i = 0, j = 0;
  const double prior = 3.0;
  int M = psis_gpd_grid_length(N);
  double xstar = x[(int) floor(N / 4.0 + 0.5) - 1];   // first quartile of sample
  double k = 0.0, sigma = 0.0, theta_hat = 0.0, l_norm = 0.0;
  LDOUBLE sum = 0.0, wt_sum = 0.0;

  if(!(xstar > x[0])){
    *k_out = NA_REAL;
    *sigma_out = NA_REAL;
    return;
  }

  // profile log-likelihood on the grid of theta values
  for(j = 0; j < M; j++){
    theta[j] = 1.0 / x[N - 1] + (1.0 - sqrt(M / ((j + 1) - 0.5))) / prior / xstar;
    sum = 0.0;
    for(i = 0; i < N; i++){
      sum += log1p(- theta[j] * x[i]);
    }
    k = (double) (sum / N);                                       // matrixStats::rowMeans2(log1p(-theta %o% x))
    l_theta[j] = N * (log(- theta[j] / k) - k - 1.0);
  }

  // posterior mean of theta with normalized weights exp(l_theta - logSumExp(l_theta))
  l_norm = psis_logSumExp(l_theta, M);
  wt_sum = 0.0;
  for(j = 0; j < M; j++){
    wt_sum += theta[j] * exp(l_theta[j] - l_norm);
  }
  theta_hat = (double) wt_sum;                                    // sum(theta * w_theta)

  for(i = 0; i < N; i++){
    x[i] = log1p(- theta_hat * x[i]);
  }
  k = psis_mean(x, N);                                            // mean.default(log1p(-theta_hat * x))
  sigma = - k / theta_hat;

  // adjust k based on the weakly informative prior, Gaussian centered on 0.5
  k = (k * N + 0.5 * 10) / (N + 10);

  if(ISNAN(k)){
    k = R_PosInf;
    sigma = R_NaN;
  }

  *k_out = k;
  *sigma_out = sigma;

}

// quantile function of the generalized Pareto distribution with location 0;
// posterior::qgeneralized_pareto(p, 0, sigma, k)
static double psis_qgpd(double p, double sigma, double k){
  if(ISNAN(sigma) || sigma <= 0) return R_NaN;
  double log_survival = log1p(-p);
  if(k == 0){
    return - sigma * log_survival;
  }
  return sigma * expm1(- k * log_survival) / k;
}

/*****************************************
 PSIS leave-one-out predictive density for one observation.

 Input:
   ll      : log-likelihood of the observation at S posterior draws
   S       : number of posterior draws
   L       : tail length, psis_tail_length(S)
 Workspace (caller allocated):
   lw      : length S; on exit, unnormalized smoothed log weights (as in loo::psis())
   idx     : length S (int)
   x_tail  : length L
   theta, l_theta : length psis_gpd_grid_length(L)
 Output:
   elpd    : log of the PSIS-LOO predictive density; loo::loo()$pointwise[, "elpd_loo"]
   khat    : Pareto k diagnostic (Inf if the tail could not be smoothed)

 Assumes ll is finite.
 *****************************************/
void psis_loo(double *ll, int S, int L, double *lw, int *idx, double *x_tail,
              double *theta, double *l_theta, double *elpd, double *khat){

  int i = 0, s = 0;
  double lr_max = R_NegInf, cutoff = 0.0, exp_cutoff = 0.0, tail_i = 0.0;
  double k = R_PosInf, sigma = 0.0, norm_const = 0.0, val_max = R_NegInf;
  LDOUBLE sum = 0.0;

  // log importance ratios are -ll; shift by their maximum for safe exponentiation
  for(s = 0; s < S; s++){
    if(- ll[s] > lr_max) lr_max = - ll[s];
  }
  for(s = 0; s < S; s++){
    lw[s] = - ll[s] - lr_max;
  }

  if(L >= 5){

    // order statistics: move the largest L + 1 log ratios to the end and sort
    // only those; ties are broken by index, as in the stable radix sort.int()
    for(s = 0; s < S; s++) idx[s] = s;
    auto by_value = [lw](int a, int b){
      return (lw[a] < lw[b]) || (lw[a] == lw[b] && a < b);
    };
    std::nth_element(idx, idx + (S - L - 1), idx + S, by_value);
    std::sort(idx + (S - L - 1), idx + S, by_value);

    cutoff = lw[idx[S - L - 1]];                                  // largest value smaller than tail values

    if(fabs(lw[idx[S - 1]] - lw[idx[S - L]]) >= DBL_EPSILON / 100){

      // exp(tail) - exp(cutoff), computed stably as -exp(x) * expm1(cutoff - x)
      for(i = 0; i < L; i++){
        tail_i = lw[idx[S - L + i]];
        x_tail[i] = (tail_i == cutoff) ? 0.0 : - exp(tail_i) * expm1(cutoff - tail_i);
      }

      psis_gpdfit(x_tail, L, theta, l_theta, &k, &sigma);
      if(ISNAN(k)) k = R_PosInf;

      // replace the tail by the expected order statistics of the fitted distribution
      if(R_FINITE(k)){
        exp_cutoff = exp(cutoff);
        for(i = 0; i < L; i++){
          lw[idx[S - L + i]] = log(psis_qgpd(((i + 1) - 0.5) / L, sigma, k) + exp_cutoff);
        }
      }

    }

  }

  // truncate at the maximum raw log ratio (0 after the shift), then shift back
  for(s = 0; s < S; s++){
    if(lw[s] > 0) lw[s] = 0.0;
    lw[s] += lr_max;
  }

  // elpd = log(sum(exp(ll + lw - logSumExp(lw))))
  norm_const = psis_logSumExp(lw, S);
  for(s = 0; s < S; s++){
    tail_i = ll[s] + (lw[s] - norm_const);
    if(tail_i > val_max) val_max = tail_i;
  }
  sum = 0.0;
  for(s = 0; s < S; s++){
    sum += exp(ll[s] + (lw[s] - norm_const) - val_max);
  }

  *elpd = val_max + log((double) sum);
  *khat = k;

}

extern "C" {

  // R interface: PSIS-LOO for an S x n matrix of log-likelihood values
  SEXP R_psis(SEXP ll_r, SEXP return_weights_r){

    int i = 0, nProtect = 0;
    int S = Rf_nrows(ll_r);
    int n = Rf_ncols(ll_r);
    int return_weights = INTEGER(return_weights_r)[0];
    double *ll = REAL(ll_r);

    int L = psis_tail_length(S);
    int M = psis_gpd_grid_length(L);

    double *lw = (double *) R_alloc(S, sizeof(double));
    int *idx = (int *) R_alloc(S, sizeof(int));
    double *x_tail = (double *) R_alloc(L > 0 ? L : 1, sizeof(double));
    double *theta = (double *) R_alloc(M, sizeof(double));
    double *l_theta = (double *) R_alloc(M, sizeof(double));

    SEXP elpd_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;
    SEXP khat_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;
    SEXP lw_r = R_NilValue;
    if(return_weights){
      lw_r = PROTECT(Rf_allocMatrix(REALSXP, S, n)); nProtect++;
    }

    for(i = 0; i < n; i++){
      psis_loo(&ll[(size_t) i * S], S, L, lw, idx, x_tail, theta, l_theta,
               &REAL(elpd_r)[i], &REAL(khat_r)[i]);
      if(return_weights){
        std::copy(lw, lw + S, &REAL(lw_r)[(size_t) i * S]);
      }
    }

    int nResultListObjs = return_weights ? 3 : 2;
    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    SEXP resultName_r = PROTECT(Rf_allocVector(STRSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, elpd_r);
    SET_STRING_ELT(resultName_r, 0, Rf_mkChar("elpd_loo"));

    SET_VECTOR_ELT(result_r, 1, khat_r);
    SET_STRING_ELT(resultName_r, 1, Rf_mkChar("pareto_k"));

    if(return_weights){
      SET_VECTOR_ELT(result_r, 2, lw_r);
      SET_STRING_ELT(resultName_r, 2, Rf_mkChar("log_weights"));
    }

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

  }

}
