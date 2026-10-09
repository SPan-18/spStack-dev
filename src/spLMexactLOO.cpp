#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
#include "psis.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

// Fit of one candidate model (fixed phi, nu and deltasq) with optional leave-one-out predictive
// densities. On entry, the lower triangle of the n x n buffer cholVy holds the spatial correlation
// matrix R; the buffer is overwritten. The prior arguments are not modified.
static SEXP spLMexactLOO_fit(double *Y, double *X, int n, int p, double *cholVy,
                             std::string &betaPrior, double *betaMu, double *betaV0,
                             double sigmaSqIGa, double sigmaSqIGb, double deltasq,
                             int nSamples, int loopd, std::string &loopd_method){

  /*****************************************
   Common variables
   *****************************************/
  int i, j, s, info, nProtect = 0;
  char const *lower = "L";
  char const *nUnit = "N";
  char const *ntran = "N";
  char const *ytran = "T";
  char const *lside = "L";
  const double one = 1.0;
  const double negOne = -1.0;
  const double zero = 0.0;
  const int incOne = 1;

  int pp = p * p;
  int np = n * p;

  const char *exact_str = "exact";
  const char *psis_str = "psis";

  // working copy of the prior covariance of beta (rescaled below)
  double *betaV = NULL;
  if(betaPrior == "normal"){
    betaV = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dcopy)(&pp, betaV0, &incOne, betaV, &incOne);
  }

  /*****************************************
   Priors on (beta, sigmaSqz) for sampling
   *****************************************/
  // The sampler below is written in terms of the spatial variance
  // sigmaSqz = sigmaSq/deltasq, where sigmaSq is the measurement error
  // variance. The prior sigmaSq ~ IG(a, b), beta | sigmaSq ~ N(muBeta, sigmaSq*Vbeta)
  // is equivalent to sigmaSqz ~ IG(a, b/deltasq), beta | sigmaSqz ~ N(muBeta, sigmaSqz*deltasq*Vbeta).
  // The flat prior p(sigmaSq) proportional to 1/sigmaSq (a = b = 0) is
  // invariant to this rescaling and so is the flat prior on beta.
  sigmaSqIGb = sigmaSqIGb / deltasq;                                                     // b = b/deltasq
  if(betaPrior == "normal"){
    F77_NAME(dscal)(&pp, &deltasq, betaV, &incOne);                                      // betaV = deltasq*Vbeta
  }

  /*****************************************
   Set-up posterior sample vector/matrices etc.
   *****************************************/
  double sigmaSqIGaPost = 0, sigmaSqIGbPost = 0;
  double sse = 0;
  double dtemp = 0;

  const double delta = sqrt(deltasq);
  char const *rside = "R";


  double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);          // allocate memory for n x 1 vector

  double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);                  // allocate memory for p x 1 vector
  double *VbetaInvMuBeta = (double *) R_alloc(p, sizeof(double)); zeros(VbetaInvMuBeta, p);  // allocate memory for p x 1 vector

  double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);  // allocate VbetaInv
  double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);      // allocate memory for p x p matrix

  // exact leave-one-out predictive densities are computed in closed form from inv(Vy)
  const int loopd_exact = (loopd && loopd_method == exact_str);
  SEXP loopd_out_r = R_NilValue;                                                   // leave-one-out predictive densities
  SEXP loopd_k_r = R_NilValue;                                                     // Pareto k diagnostics (PSIS only)
  if(loopd){
    loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;
  }
  double *looQr = NULL;                                                            // Q*r = inv(Vy)*(Y-X*betahat), exact LOO only

  // diagnostics: correlations of the farthest-apart and the closest locations, and below the smallest relative
  // Cholesky pivot of the two n x n factorizations (no extra factorization)
  double diagMinCor = 0.0, diagMaxCor = 0.0, diagPivot = 0.0;
  corOffDiagRange(cholVy, n, &diagMinCor, &diagMaxCor);

  // construct marginal covariance matrix Vy = R + deltasq*I in place
  for(i = 0; i < n; i++){
    cholVy[i*n + i] += deltasq;
  }

  // chol(Vy)
  F77_NAME(dpotrf)(lower, &n, cholVy, &n, &info FCONE);
  if(info != 0){Rf_error("c++ error: Cholesky factorization of Vy = R + deltasq*I failed (info = %i).\n", info);}
  diagPivot = minRelPivot(cholVy, n, NULL, 1.0 + deltasq);                                  // diag(Vy) = 1 + deltasq

  // find cholinv(Vy)*Y
  F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);                                         // tmp_n = Y
  F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVy, &n, tmp_n, &incOne FCONE FCONE FCONE);  // tmp_n = cholinv(Vy)*Y

  if(betaPrior == "normal"){
    // find VbetaInvmuBeta
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                     // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE);                                                     // VbetaInv = chol(Vbeta)
    if(info != 0){Rf_error("c++ error: prior covariance of beta is not positive definite.\n");}
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE);                                                     // VbetaInv = chol2inv(Vbeta)
    if(info != 0){Rf_error("c++ error: inversion of the prior covariance of beta failed.\n");}
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betaMu, &incOne, &zero, VbetaInvMuBeta, &incOne FCONE);       // VbetaInvMuBeta = VbetaInv*muBeta
  }else{
    // flat prior on beta: VbetaInv = 0 and VbetaInv*muBeta = 0 (already
    // zero-initialized); since p(beta) does not carry the factor
    // (sigmaSq)^(-p/2) of the conjugate prior, the posterior inverse-gamma
    // shape is aIG + (n - p)/2
    sigmaSqIGa -= 0.5 * p;
  }

  //  find XtVyInvY
  double *tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np, np);                            // allocate temporary memory for n x p matrix
  F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);                                                          // tmp_np = X
  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, cholVy, &n, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = cholinv(Vy)*X
  F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);             // tmp_p1 = t(X)*VyInv*Y

  // find betahat = inv(XtVyInvX + VbetaInv)(XtVyInvY + VbetaInvmuBeta)
  F77_NAME(daxpy)(&p, &one, VbetaInvMuBeta, &incOne, tmp_p1, &incOne);                                        // tmp_p1 = XtVyInvY + VbetaInvmuBeta
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);     // tmp_pp = t(X)*VyInv*X

  F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, tmp_pp, &incOne);                                             // tmp_pp = t(X)*VyInv*X + VbetaInv
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);                                                       // tmp_pp = chol(XtVyInvX + VbetaInv)
  if(info != 0){
    R_chk_free(tmp_np);
    Rf_error("c++ error: Cholesky factorization of t(X)*inv(Vy)*X + inv(Vbeta) failed; check X for collinear columns.\n");
  }
  F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                    // tmp_p1 = cholinv(XtVyInvX + VbetaInv)*tmp_p1

  // find sse = t(Y-X*betahat)*VyInv*(Y-X*betahat) + t(betahat-muBeta)*VbetaInv*(betahat-muBeta);
  // equal to t(Y)*VyInv*Y + t(muBeta)*VbetaInv*muBeta - t(m)*M*m, but each term is a sum of
  // squares, hence non-negative and free of catastrophic cancellation
  double *betahat = (double *) R_chk_calloc(p, sizeof(double)); zeros(betahat, p);                           // allocate temporary memory for p x 1 vector
  F77_NAME(dcopy)(&p, tmp_p1, &incOne, betahat, &incOne);                                                     // betahat = cholinv(XtVyInvX + VbetaInv)*tmp_p1
  F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, betahat, &incOne FCONE FCONE FCONE);                   // betahat = inv(XtVyInvX + VbetaInv)*(XtVyInvY + VbetaInvmuBeta)
  F77_NAME(dgemv)(ntran, &n, &p, &negOne, tmp_np, &n, betahat, &incOne, &one, tmp_n, &incOne FCONE);          // tmp_n = cholinv(Vy)*(Y-X*betahat)
  sse = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne), 2);                                                          // sse = t(Y-X*betahat)*VyInv*(Y-X*betahat)

  if(loopd_exact){
    // ingredients of the closed-form leave-one-out predictive densities, see below
    looQr = (double *) R_chk_calloc(n, sizeof(double)); zeros(looQr, n);                                      // allocate temporary memory for n x 1 vector
    F77_NAME(dcopy)(&n, tmp_n, &incOne, looQr, &incOne);                                                      // looQr = cholinv(Vy)*(Y-X*betahat)
    F77_NAME(dtrsv)(lower, ytran, nUnit, &n, cholVy, &n, looQr, &incOne FCONE FCONE FCONE);                   // looQr = inv(Vy)*(Y-X*betahat) = Q*r
    F77_NAME(dtrsm)(lside, lower, ytran, nUnit, &n, &p, &one, cholVy, &n, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = inv(Vy)*X
    F77_NAME(dtrsm)(rside, lower, ytran, nUnit, &n, &p, &one, tmp_pp, &p, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = inv(Vy)*X*t(cholinv(XtVyInvX + VbetaInv))
  }

  if(betaPrior == "normal"){
    F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, betahat, &incOne);                                          // betahat = betahat - muBeta
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betahat, &incOne, &zero, tmp_n, &incOne FCONE);            // tmp_n[1:p] = VbetaInv*(betahat - muBeta)
    sse += F77_CALL(ddot)(&p, betahat, &incOne, tmp_n, &incOne);                                              // sse = sse + t(betahat-muBeta)*VbetaInv*(betahat-muBeta)
  }

  // deallocate betahat
  R_chk_free(betahat);

  if(!(sigmaSqIGb + 0.5 * sse > 0.0)){
    R_chk_free(tmp_np);
    if(loopd_exact){ R_chk_free(looQr); }
    Rf_error("c++ error: posterior scale of sigmaSq is not positive (sse = %g).\n", sse);
  }

  // inv(Vy) from chol(Vy), in place (lower triangle only)
  F77_NAME(dpotri)(lower, &n, cholVy, &n, &info FCONE);                                                       // cholVy = inv(Vy)
  if(info != 0){
    R_chk_free(tmp_np);
    if(loopd_exact){ R_chk_free(looQr); }
    Rf_error("c++ error: inversion of Vy failed (info = %i).\n", info);
  }

  /*****************************************
   Exact leave-one-out predictive densities
   *****************************************/
  // Let K = Vy + X*Vbeta*t(X) (backend units), Q = inv(K) = inv(Vy) - inv(Vy)*X*B*t(X)*inv(Vy) with
  // B = inv(t(X)*inv(Vy)*X + inv(Vbeta)) (inv(Vbeta) = 0 under a flat prior on beta), and r = Y - X*muBeta.
  // Marginally, Y is multivariate-t; by the partitioned-inverse identities, the predictive density of
  // Y[i] given Y[-i] (the exact LOO-PD expression of the paper, with every inv(V_{y_{-i}}) term written through Q) is
  //   t_{2a_i}(Y[i]; Y[i] - (Qr)_i/Q_ii, (b_i/a_i)/Q_ii),
  //   a_i = a + (n - 1)/2  (a already reduced by p/2 under a flat prior on beta),
  //   b_i = b + (t(r)*Q*r - (Qr)_i^2/Q_ii)/2,
  // with Q*r = inv(Vy)*(Y - X*betahat), t(r)*Q*r = sse, and Q_ii = inv(Vy)_ii - ||row i of inv(Vy)*X*t(cholinv(B^-1))||^2.
  if(loopd_exact){

    double a_i = sigmaSqIGa + 0.5 * (n - 1);
    double Q_ii = 0.0, b_i = 0.0, scale = 0.0;

    for(i = 0; i < n; i++){
      Q_ii = cholVy[i*n + i];                                                                                 // Q_ii = inv(Vy)_ii
      for(j = 0; j < p; j++){
        Q_ii -= tmp_np[j*n + i] * tmp_np[j*n + i];                                                            // Q_ii = inv(Vy)_ii - t(w_i)*B*w_i
      }
      dtemp = looQr[i] / Q_ii;                                                                                // dtemp = (Qr)_i/Q_ii = Y[i] - location_i
      b_i = sse - looQr[i] * dtemp;                                                                           // b_i = t(r)*Q*r - (Qr)_i^2/Q_ii >= 0
      b_i = sigmaSqIGb + 0.5 * fmax2(b_i, 0.0);
      scale = sqrt((b_i / a_i) / Q_ii);
      REAL(loopd_out_r)[i] = Rf_dt(dtemp / scale, 2.0 * a_i, 1) - log(scale);
    }

    R_chk_free(looQr);

  }

  // deallocate tmp_np
  R_chk_free(tmp_np);

  /*****************************************
   Set-up for sampling spatial random effects
   *****************************************/
  // the posterior covariance of z is sigmaSq*deltasq*inv(Vy)*R; since Vy = R + deltasq*I,
  // inv(Vy)*R = I - deltasq*inv(Vy), formed in place (lower triangle only)
  const double negDeltasq = -1.0 * deltasq;
  for(j = 0; j < n; j++){
    i = n - j;
    F77_NAME(dscal)(&i, &negDeltasq, &cholVy[j*n + j], &incOne);                                            // cholVy[j:n, j] = -deltasq*inv(Vy)[j:n, j]
    cholVy[j*n + j] += 1.0;                                                                                   // cholVy = I - deltasq*inv(Vy) = inv(Vy)*R
  }
  const int nPlusOne = n + 1;
  double *diagM = (double *) R_alloc(n, sizeof(double));
  F77_NAME(dcopy)(&n, cholVy, &nPlusOne, diagM, &incOne);                                                     // diag(inv(Vy)*R)
  F77_NAME(dpotrf)(lower, &n, cholVy, &n, &info FCONE);                                                       // cholVy = chol(inv(Vy)*R)
  if(info != 0){Rf_error("c++ error: Cholesky factorization of the posterior covariance of z failed (info = %i); the spatial correlation matrix is numerically singular.\n", info);}
  diagPivot = fmin2(diagPivot, minRelPivot(cholVy, n, diagM, 1.0));

  // posterior parameters of sigmaSq
  sigmaSqIGaPost += sigmaSqIGa;
  sigmaSqIGaPost += 0.5 * n;

  sigmaSqIGbPost += sigmaSqIGb;
  sigmaSqIGbPost += 0.5 * sse;

  // posterior samples of sigma-sq and beta
  SEXP samples_sigmaSq_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
  SEXP samples_sigmaSqz_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
  SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
  SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

  // Composition sampling: for each s, draw sigmaSqz_s, then beta_s, then
  //   z_s = L*(xi_s + t(L)*(Y - X*beta_s)),  xi_s ~ N(0, deltasq*sigmaSqz_s*I),  L = chol(inv(Vy)*R),
  // so that z_s ~ N(inv(Vy)*R*(Y - X*beta_s), sigmaSqz_s*deltasq*inv(Vy)*R). The random variates are
  // drawn in the loop (xi_s stored in column s of samples_z); the linear algebra is then done for all
  // draws at once: Z = L*(Xi + t(L)*Y*1' - t(L)*X*B), with level-3 BLAS and no extra n x nSamples storage.
  double sigmaSqz = 0;
  double *pointer_beta = REAL(samples_beta_r);
  double *pointer_z = REAL(samples_z_r);
  double *beta_s = NULL, *z_s = NULL;

  GetRNGstate();

  for(s = 0; s < nSamples; s++){
    // sample sigmaSqz (spatial variance) from its marginal posterior
    dtemp = 1.0 / sigmaSqIGbPost;
    dtemp = rgamma(sigmaSqIGaPost, dtemp);
    sigmaSqz = 1.0 / dtemp;
    REAL(samples_sigmaSqz_r)[s] = sigmaSqz;
    REAL(samples_sigmaSq_r)[s] = deltasq * sigmaSqz;                                       // sigmaSq = deltasq*sigmaSqz

    // sample fixed effects by composition sampling
    beta_s = &pointer_beta[(R_xlen_t) s * p];
    dtemp = sqrt(sigmaSqz);
    for(j = 0; j < p; j++){
      beta_s[j] = rnorm(tmp_p1[j], dtemp);                                                 // beta ~ N(tmp_p1, sigmaSq*I)
    }
    F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, beta_s, &incOne FCONE FCONE FCONE); // beta = t(cholinv(tmp_pp))*beta

    // random part of the spatial effects
    z_s = &pointer_z[(R_xlen_t) s * n];
    dtemp = dtemp * delta;                                                                 // dtemp = sqrt(deltasq*sigmaSqz)
    for(i = 0; i < n; i++){
      z_s[i] = rnorm(0.0, dtemp);                                                          // xi_s ~ N(0, deltasq*sigmaSqz*I)
    }

  }

  PutRNGstate();

  // spatial effects for all draws: Z = L*(Xi + t(L)*Y*1' - t(L)*X*B), L = chol(inv(Vy)*R) (lower triangle)
  double *tmp_np2 = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np2, np);                         // allocate temporary memory for n x p matrix
  F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);                                                           // tmp_n = Y
  F77_NAME(dtrmv)(lower, ytran, nUnit, &n, cholVy, &n, tmp_n, &incOne FCONE FCONE FCONE);                    // tmp_n = t(L)*Y
  F77_NAME(dcopy)(&np, X, &incOne, tmp_np2, &incOne);                                                        // tmp_np2 = X
  F77_NAME(dtrmm)(lside, lower, ytran, nUnit, &n, &p, &one, cholVy, &n, tmp_np2, &n FCONE FCONE FCONE FCONE); // tmp_np2 = t(L)*X
  for(s = 0; s < nSamples; s++){
    F77_NAME(daxpy)(&n, &one, tmp_n, &incOne, &pointer_z[(R_xlen_t) s * n], &incOne);                       // Z[, s] = xi_s + t(L)*Y
  }
  F77_NAME(dgemm)(ntran, ntran, &n, &nSamples, &p, &negOne, tmp_np2, &n, pointer_beta, &p, &one, pointer_z, &n FCONE FCONE); // Z = Xi + t(L)*Y*1' - t(L)*X*B
  F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &n, &nSamples, &one, cholVy, &n, pointer_z, &n FCONE FCONE FCONE FCONE);    // Z = L*Z
  R_chk_free(tmp_np2);

  // make return object
  SEXP result_r, resultName_r;

  // If loopd is TRUE, set-up Leave-one-out predictive density calculation
  if(loopd){

    if(loopd_method == psis_str){

      int loo_index = 0;
      double theta_i = 0.0;

      // PSIS workspace, O(nSamples), allocated once for all observations
      int psis_L = psis_tail_length(nSamples);
      int psis_M = psis_gpd_grid_length(psis_L);

      double *X_i = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_i, p);
      double *ll_i = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(ll_i, nSamples);
      double *lw_i = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(lw_i, nSamples);
      int *idx_i = (int *) R_chk_calloc(nSamples, sizeof(int)); zeros(idx_i, nSamples);
      double *xtail_i = (double *) R_chk_calloc(psis_L, sizeof(double)); zeros(xtail_i, psis_L);
      double *theta_gpd = (double *) R_chk_calloc(psis_M, sizeof(double)); zeros(theta_gpd, psis_M);
      double *ltheta_gpd = (double *) R_chk_calloc(psis_M, sizeof(double)); zeros(ltheta_gpd, psis_M);

      loopd_k_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

      double *pointer_sigmaSq = REAL(samples_sigmaSq_r);

      for(loo_index = 0; loo_index < n; loo_index++){

        copyMatrixRowToVec(X, n, p, X_i, loo_index);                                          // X_i = X[i,1:p]

        // log-likelihood of the i-th observation at each posterior draw
        for(s = 0; s < nSamples; s++){
          theta_i = F77_CALL(ddot)(&p, X_i, &incOne, &pointer_beta[(R_xlen_t) s * p], &incOne);            // theta_i = X_i * beta_s
          theta_i += pointer_z[(R_xlen_t) s * n + loo_index];                                              // theta_i = X_i*beta_s + zi_s
          ll_i[s] = Rf_dnorm4(Y[loo_index], theta_i, sqrt(pointer_sigmaSq[s]), 1);           // sigmaSq_s is the measurement error variance
        }

        psis_loo(ll_i, nSamples, psis_L, lw_i, idx_i, xtail_i, theta_gpd, ltheta_gpd,
                 &REAL(loopd_out_r)[loo_index], &REAL(loopd_k_r)[loo_index]);

      }

      R_chk_free(X_i);
      R_chk_free(ll_i);
      R_chk_free(lw_i);
      R_chk_free(idx_i);
      R_chk_free(xtail_i);
      R_chk_free(theta_gpd);
      R_chk_free(ltheta_gpd);

    }

    // make return object for posterior samples and leave-one-out predictive densities
    int nResultListObjs = (loopd_method == psis_str) ? 6 : 5;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

    // samples of sigma-sq
    SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

    // samples of sigma-sq-z
    SET_VECTOR_ELT(result_r, 2, samples_sigmaSqz_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("sigmaSq.z"));

    // samples of z
    SET_VECTOR_ELT(result_r, 3, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z"));

    // leave-one-out predictive densities
    SET_VECTOR_ELT(result_r, 4, loopd_out_r);
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("loopd"));

    // Pareto k diagnostics of PSIS
    if(loopd_method == psis_str){
      SET_VECTOR_ELT(result_r, 5, loopd_k_r);
      SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("loopd.pareto_k"));
    }

    Rf_namesgets(result_r, resultName_r);

  }else{

    // make return object for posterior samples of sigma-sq, beta and z
    int nResultListObjs = 4;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

    // samples of sigma-sq
    SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

    // samples of sigma-sq-z
    SET_VECTOR_ELT(result_r, 2, samples_sigmaSqz_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("sigmaSq.z"));

    // samples of z
    SET_VECTOR_ELT(result_r, 3, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z"));

    Rf_namesgets(result_r, resultName_r);

  }

  result_r = PROTECT(appendDiagnostics(result_r, diagPivot, diagMinCor, diagMaxCor)); nProtect++;

  UNPROTECT(nProtect);

  return result_r;
}

extern "C" {

  SEXP spLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coords_r,
                    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                    SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                    SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                    SEXP verbose_r){

    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int n = INTEGER(n_r)[0];
    int nn = n * n;

    // Set-up coordinates (n x 2) and spatial correlation function
    double *coords = REAL(coords_r);
    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r, 0));
    double *betaMu = NULL;
    double *betaV = NULL;

    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);

      betaV = (double *) R_alloc(pp, sizeof(double));
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaV, &incOne);
    }

    double sigmaSqIGa = REAL(sigmaSqIG_r)[0];
    double sigmaSqIGb = REAL(sigmaSqIG_r)[1];

    double deltasq = REAL(deltasq_r)[0];
    double phi = REAL(phi_r)[0];

    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    // Leave-one-out predictive density details
    int loopd = INTEGER(loopd_r)[0];
    std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];

    // print set-up if verbose TRUE
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept).\n\n", p);
      Rprintf("Using the %s spatial correlation function.\n\n", corfn.c_str());

      Rprintf("Priors:\n");

      if(betaPrior == "flat"){
        Rprintf("\tbeta flat.\n");
      }else{
        Rprintf("\tbeta: Gaussian\n");
        Rprintf("\tmu:"); printVec(betaMu, p);
        Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
        Rprintf("\n");
      }

      // prior on the measurement error variance sigma.sq; b = 0 corresponds
      // to the flat prior p(sigma.sq) proportional to 1/sigma.sq
      if(sigmaSqIGb == 0.0){
        Rprintf("\tsigma.sq: flat, proportional to 1/sigma.sq.\n\n");
      }else{
        Rprintf("\tsigma.sq: Inverse-Gamma\n\tshape = %.2f, scale = %.2f.\n\n",
                sigmaSqIGa, sigmaSqIGb);
      }

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.2f, and, nu = %.2f.\n", phi, nu);
      }else{
        Rprintf("\tphi = %.2f.\n", phi);
      }
      Rprintf("Noise-to-spatial variance ratio = %.2f.\n\n", deltasq);

      Rprintf("Number of posterior samples = %i.\n\n", nSamples);

      if(loopd){
        Rprintf("LOO-PD calculation method = %s.\n", loopd_method.c_str());
      }

      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Spatial correlation matrix and model fit
     *****************************************/
    double *cholVy = (double *) R_alloc(nn, sizeof(double)); zeros(cholVy, nn);      // the only n x n matrix: R -> Vy -> chol(Vy) -> inv(Vy) -> chol(inv(Vy)*R)
    double thetasp[2] = {phi, nu};                                                   // spatial process parameters
    spCorFull2(n, 2, coords, thetasp, corfn, cholVy);

    return spLMexactLOO_fit(Y, X, n, p, cholVy, betaPrior, betaMu, betaV,
                            sigmaSqIGa, sigmaSqIGb, deltasq, nSamples, loopd, loopd_method);

  }


  // Fits of the candidate models sharing (phi, nu) for a vector of deltasq values. The correlation
  // matrix R is built once and kept in packed lower-triangular storage (n(n+1)/2 doubles); for each
  // deltasq its lower triangle is restored into the n x n working buffer (all routines read only the
  // lower triangle). Returns a list with one element per deltasq, each as returned by spLMexactLOO.
  SEXP spLMexactLOOgrid(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coords_r,
                        SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                        SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                        SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r){

    int j, k, len;
    size_t offset = 0;
    const int incOne = 1;

    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int n = INTEGER(n_r)[0];
    int nn = n * n;

    double *coords = REAL(coords_r);
    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r, 0));
    double *betaMu = NULL;
    double *betaV = NULL;
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
      betaV = (double *) R_alloc(pp, sizeof(double));
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaV, &incOne);
    }
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0];
    double sigmaSqIGb = REAL(sigmaSqIG_r)[1];

    double phi = REAL(phi_r)[0];
    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    int nModels = Rf_length(deltasq_r);
    double *deltasq = REAL(deltasq_r);

    int nSamples = INTEGER(nSamples_r)[0];
    int loopd = INTEGER(loopd_r)[0];
    std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));

    // correlation matrix R, built once
    double *cholVy = (double *) R_alloc(nn, sizeof(double)); zeros(cholVy, nn);      // n x n working buffer
    double thetasp[2] = {phi, nu};                                                   // spatial process parameters
    spCorFull2(n, 2, coords, thetasp, corfn, cholVy);

    size_t nPacked = (size_t) n * (n + 1) / 2;
    double *Rpacked = (double *) R_alloc(nPacked, sizeof(double));                   // lower triangle of R, packed by columns
    offset = 0;
    for(j = 0; j < n; j++){
      len = n - j;
      F77_NAME(dcopy)(&len, &cholVy[j*n + j], &incOne, &Rpacked[offset], &incOne);  // Rpacked <- R[j:n, j]
      offset += len;
    }

    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, nModels));

    for(k = 0; k < nModels; k++){

      const void *vmax = vmaxget();                                                  // release the per-model R_alloc memory below

      offset = 0;
      for(j = 0; j < n; j++){
        len = n - j;
        F77_NAME(dcopy)(&len, &Rpacked[offset], &incOne, &cholVy[j*n + j], &incOne);  // cholVy[j:n, j] <- R[j:n, j]
        offset += len;
      }

      SET_VECTOR_ELT(result_r, k, spLMexactLOO_fit(Y, X, n, p, cholVy, betaPrior, betaMu, betaV,
                                                   sigmaSqIGa, sigmaSqIGb, deltasq[k], nSamples,
                                                   loopd, loopd_method));
      vmaxset(vmax);

      R_CheckUserInterrupt();

    }

    UNPROTECT(1);

    return result_r;

  }

}
