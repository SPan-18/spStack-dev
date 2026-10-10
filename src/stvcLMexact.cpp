#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include <cstring>
#include "util.h"
#include "psis.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

// Bayesian linear model with spatially-temporally varying coefficients (Gaussian response),
//   y = X*beta + sum_q D_q*z_q + eps,  eps ~ N(0, sigmaSq*I),  D_q = diag(XTilde[, q]),
//   z_q | sigmaSq ~ N(0, (sigmaSq/deltasq_q)*R_q) independent,
//   beta | sigmaSq ~ N(muBeta, sigmaSq*Vbeta) (or flat),  sigmaSq ~ IG(a, b) (or p(sigmaSq) ~ 1/sigmaSq),
// with the noise-to-spatial variance ratios deltasq_q and the correlation matrices R_q = R(phi_s_q, phi_t_q)
// fixed. With Rcal = blockdiag(R_q/deltasq_q) and XTildeBlk = [D_1, ..., D_r], the marginal covariance of y given
// (beta, sigmaSq) is sigmaSq*Vy, Vy = I + XTildeBlk*Rcal*t(XTildeBlk) = I + sum_q D_q*R_q*D_q/deltasq_q, and the joint
// posterior is sampled exactly by composition:
//   sigmaSq | y ~ IG(a + n/2, b + sse/2),  beta | sigmaSq, y ~ N(betahat, sigmaSq*B),
//   z | beta, sigmaSq, y ~ N(Rcal*t(XTildeBlk)*inv(Vy)*(y - X*beta), sigmaSq*(Rcal - Rcal*t(XTildeBlk)*inv(Vy)*XTildeBlk*Rcal)).
// z is drawn with the update of Bhattacharya, Chakraborty and Mallick (2016) (Matheron's rule): with
// u ~ N(0, sigmaSq*Rcal) and d ~ N(0, sigmaSq*I),
//   z = u + Rcal*t(XTildeBlk)*inv(Vy)*(y - X*beta - XTildeBlk*u - d)
// has exactly the distribution above. This needs chol(Vy) and chol(R_q) (n x n), never the nr x nr posterior
// covariance of z, and all draws are processed with level-3 BLAS.
// Storage: for each distinct correlation matrix, one n x n buffer holds chol(R_q) in its lower triangle and R_q in
// its strict upper triangle (diag(R_q) = 1); one n x n buffer holds the upper Cholesky factor U of Vy = t(U)*U,
// replaced by inv(U) if exact leave-one-out predictive densities are required (they need diag(inv(Vy))).

// Upper Cholesky factor U of Vy, or its inverse (inv = 1)
struct VyFactor {
  double *U;
  int n;
  int inv;
};

// B <- inv(t(U))*B (B is n x ncol with leading dimension ldb)
static void solveUt(VyFactor &F, int ncol, double *B, int ldb){
  const double one = 1.0;
  if(F.inv){
    F77_NAME(dtrmm)("L", "U", "T", "N", &F.n, &ncol, &one, F.U, &F.n, B, &ldb FCONE FCONE FCONE FCONE);
  }else{
    F77_NAME(dtrsm)("L", "U", "T", "N", &F.n, &ncol, &one, F.U, &F.n, B, &ldb FCONE FCONE FCONE FCONE);
  }
}

// B <- inv(U)*B
static void solveU(VyFactor &F, int ncol, double *B, int ldb){
  const double one = 1.0;
  if(F.inv){
    F77_NAME(dtrmm)("L", "U", "N", "N", &F.n, &ncol, &one, F.U, &F.n, B, &ldb FCONE FCONE FCONE FCONE);
  }else{
    F77_NAME(dtrsm)("L", "U", "N", "N", &F.n, &ncol, &one, F.U, &F.n, B, &ldb FCONE FCONE FCONE FCONE);
  }
}

// Builds the nR distinct spatial-temporal correlation matrices R_k in Rbuf (nR blocks of n x n) and factorizes them
// in place: on exit, the lower triangle of block k holds chol(R_k) and its strict upper triangle R_k. Also records
// the off-diagonal correlation range and the smallest relative Cholesky pivot of each R_k (diagnostics).
static void stvcLM_corChol(int n, int nR, double *coords_sp, double *coords_tm, double *phi_s, double *phi_t,
                           std::string &corfn, double *Rbuf, double *pivR, double *minCor, double *maxCor){

  int k, info = 0;
  const size_t nn = (size_t) n * n;
  double thetaspt[2] = {0.0, 0.0};
  double *Rk = NULL;

  for(k = 0; k < nR; k++){
    Rk = Rbuf + nn * k;
    thetaspt[0] = phi_s[k];
    thetaspt[1] = phi_t[k];
    sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, Rk);
    corOffDiagRange(Rk, n, &minCor[k], &maxCor[k]);
    F77_NAME(dpotrf)("L", &n, Rk, &n, &info FCONE);
    if(info != 0){
      if(nR > 1){
        Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix of process %i failed (info = %i); check for nearly coincident space-time locations or very small decay parameters.\n", k + 1, info);
      }
      Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i); check for nearly coincident space-time locations or very small decay parameters.\n", info);
    }
    pivR[k] = minRelPivot(Rk, n, NULL, 1.0);
  }

}

// Fit of one candidate model (fixed correlation matrices and deltasq[0..r-1]) with optional leave-one-out predictive
// densities ('exact' in closed form, or 'psis'). Rbuf is as returned by stvcLM_corChol and is only read, so the
// candidate models sharing (phi_s, phi_t) use the same Rbuf. nR = 1 (shared correlation matrix and deltasq) or r.
static SEXP stvcLMexact_fit(double *Y, double *X, double *XTilde, int n, int p, int r, int nR,
                            double *Rbuf, double *pivR, double *minCor, double *maxCor,
                            std::string &betaPrior, double *betaMu, double *betaV,
                            double sigmaSqIGa, double sigmaSqIGb, double *deltasq,
                            int nSamples, int loopd, std::string &loopd_method){

  int i, j, k, q, s, b, info, nProtect = 0;
  const char *lower = "L";
  const char *ntran = "N";
  const char *ytran = "T";
  const char *nUnit = "N";
  const char *lside = "L";
  const char *rside = "R";
  const double one = 1.0;
  const double negOne = -1.0;
  const double zero = 0.0;
  const int incOne = 1;

  const size_t nn = (size_t) n * n;
  const int nr = n * r;
  const int pp = p * p;
  const int np = n * p;

  const int loopd_exact = (loopd && loopd_method == "exact");
  const int loopd_psis = (loopd && loopd_method == "psis");

  double dtemp = 0.0;
  double *Rk = NULL, *xq = NULL;

  /*****************************************
   Vy = I + sum_q D_q*R_q*D_q/deltasq_q and its Cholesky factor
   *****************************************/
  // upper triangle, column by column (R_q is read from the strict upper triangle of its buffer, diag(R_q) = 1)
  double *U = (double *) R_alloc(nn, sizeof(double));
  std::memset(U, 0, sizeof(double) * nn);
  double *dVy = (double *) R_alloc(n, sizeof(double));
  for(q = 0; q < r; q++){
    Rk = Rbuf + nn * (nR == 1 ? 0 : q);
    xq = XTilde + (size_t) n * q;
    const double w = 1.0 / deltasq[q];
    for(j = 0; j < n; j++){
      const double xj = xq[j] * w;
      double *Uc = U + (size_t) j * n;
      const double *Rc = Rk + (size_t) j * n;
      for(i = 0; i < j; i++){
        Uc[i] += xq[i] * xj * Rc[i];
      }
      Uc[j] += xq[j] * xj;
    }
  }
  for(i = 0; i < n; i++){
    U[(size_t) i * n + i] += 1.0;
    dVy[i] = U[(size_t) i * n + i];
  }

  F77_NAME(dpotrf)("U", &n, U, &n, &info FCONE);
  if(info != 0){Rf_error("c++ error: Cholesky factorization of Vy = I + sum_q D_q*R_q*D_q/deltasq_q failed (info = %i).\n", info);}
  const double pivVy = minRelPivot(U, n, dVy, 1.0);

  VyFactor F = {U, n, 0};

  // exact leave-one-out predictive densities need diag(inv(Vy)) = row sums of squares of inv(U); U is replaced by
  // inv(U) and all the solves below become triangular multiplications (same cost)
  double *dinvVy = NULL;
  if(loopd_exact){
    F77_NAME(dtrtri)("U", nUnit, &n, U, &n, &info FCONE FCONE);
    if(info != 0){Rf_error("c++ error: inversion of the Cholesky factor of Vy failed (info = %i).\n", info);}
    F.inv = 1;
    dinvVy = (double *) R_alloc(n, sizeof(double)); zeros(dinvVy, n);
    for(j = 0; j < n; j++){
      const double *Uc = U + (size_t) j * n;
      for(i = 0; i <= j; i++){
        dinvVy[i] += Uc[i] * Uc[i];
      }
    }
  }

  /*****************************************
   Marginal posterior of (beta, sigmaSq)
   *****************************************/
  double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);
  double *VbetaInvMuBeta = (double *) R_alloc(p, sizeof(double)); zeros(VbetaInvMuBeta, p);
  if(betaPrior == "normal"){
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: prior covariance of beta is not positive definite.\n");}
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: inversion of the prior covariance of beta failed.\n");}
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betaMu, &incOne, &zero, VbetaInvMuBeta, &incOne FCONE);
  }else{
    // flat prior on beta: p(beta) does not carry the factor sigmaSq^(-p/2) of the conjugate prior
    sigmaSqIGa -= 0.5 * p;
  }

  double *tmp_n = (double *) R_alloc(n, sizeof(double));
  double *tmp_np = (double *) R_alloc(np, sizeof(double));
  double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);
  double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);
  double *betahat = (double *) R_alloc(p, sizeof(double)); zeros(betahat, p);

  F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);
  solveUt(F, 1, tmp_n, n);                                                                                     // tmp_n = inv(t(U))*Y
  F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);
  solveUt(F, p, tmp_np, n);                                                                                    // tmp_np = inv(t(U))*X
  F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);              // t(X)*inv(Vy)*Y
  F77_NAME(daxpy)(&p, &one, VbetaInvMuBeta, &incOne, tmp_p1, &incOne);                                         // + inv(Vbeta)*muBeta
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);      // t(X)*inv(Vy)*X
  F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, tmp_pp, &incOne);                                              // + inv(Vbeta) = inv(B)
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);                                                        // tmp_pp = chol(inv(B))
  if(info != 0){Rf_error("c++ error: Cholesky factorization of t(X)*inv(Vy)*X + inv(Vbeta) failed; check X for collinear columns.\n");}
  F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                     // tmp_p1 = cholinv(inv(B))*b
  F77_NAME(dcopy)(&p, tmp_p1, &incOne, betahat, &incOne);
  F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, betahat, &incOne FCONE FCONE FCONE);                    // betahat = B*b
  F77_NAME(dgemv)(ntran, &n, &p, &negOne, tmp_np, &n, betahat, &incOne, &one, tmp_n, &incOne FCONE);           // tmp_n = inv(t(U))*(Y - X*betahat)
  double sse = F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n, &incOne);                                             // t(Y - X*betahat)*inv(Vy)*(Y - X*betahat)

  double *looQr = NULL;
  if(loopd_exact){
    looQr = (double *) R_alloc(n, sizeof(double));
    F77_NAME(dcopy)(&n, tmp_n, &incOne, looQr, &incOne);
    solveU(F, 1, looQr, n);                                                                                    // looQr = inv(Vy)*(Y - X*betahat)
    solveU(F, p, tmp_np, n);                                                                                   // tmp_np = inv(Vy)*X
    F77_NAME(dtrsm)(rside, lower, ytran, nUnit, &n, &p, &one, tmp_pp, &p, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = inv(Vy)*X*t(cholinv(inv(B)))
  }

  if(betaPrior == "normal"){
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betahat, &incOne, tmp_p2, &incOne);
    F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, tmp_p2, &incOne);                                            // betahat - muBeta
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, tmp_p2, &incOne, &zero, tmp_n, &incOne FCONE);              // tmp_n[1:p] = inv(Vbeta)*(betahat - muBeta)
    sse += F77_NAME(ddot)(&p, tmp_p2, &incOne, tmp_n, &incOne);
  }

  if(!(sigmaSqIGb + 0.5 * sse > 0.0)){
    Rf_error("c++ error: posterior scale of sigmaSq is not positive (sse = %g).\n", sse);
  }
  const double sigmaSqIGaPost = sigmaSqIGa + 0.5 * n;
  const double sigmaSqIGbPost = sigmaSqIGb + 0.5 * sse;

  /*****************************************
   Exact leave-one-out predictive densities
   *****************************************/
  // As for the spatial linear model: Y is marginally multivariate-t with scale matrix K = Vy + X*Vbeta*t(X); with
  // Q = inv(K) = inv(Vy) - inv(Vy)*X*B*t(X)*inv(Vy) and r = Y - X*muBeta, the predictive density of Y[i] given Y[-i] is
  //   t_{2a_i}(Y[i]; Y[i] - (Qr)_i/Q_ii, (b_i/a_i)/Q_ii),  a_i = a + (n - 1)/2,  b_i = b + (t(r)*Q*r - (Qr)_i^2/Q_ii)/2,
  // with Q*r = inv(Vy)*(Y - X*betahat), t(r)*Q*r = sse and Q_ii = inv(Vy)_ii - ||row i of inv(Vy)*X*t(cholinv(inv(B)))||^2.
  SEXP loopd_out_r = R_NilValue, loopd_k_r = R_NilValue;
  if(loopd){
    loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;
  }
  if(loopd_exact){
    const double a_i = sigmaSqIGa + 0.5 * (n - 1);
    double Q_ii = 0.0, b_i = 0.0, scale = 0.0;
    for(i = 0; i < n; i++){
      Q_ii = dinvVy[i];
      for(j = 0; j < p; j++){
        Q_ii -= tmp_np[(size_t) j * n + i] * tmp_np[(size_t) j * n + i];
      }
      dtemp = looQr[i] / Q_ii;
      b_i = sigmaSqIGb + 0.5 * fmax2(sse - looQr[i] * dtemp, 0.0);
      scale = sqrt((b_i / a_i) / Q_ii);
      REAL(loopd_out_r)[i] = Rf_dt(dtemp / scale, 2.0 * a_i, 1) - log(scale);
    }
  }

  /*****************************************
   Posterior samples
   *****************************************/
  SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
  SEXP samples_sigmaSq_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
  SEXP samples_sigmaSqz_r = R_NilValue;
  if(nR == 1){
    samples_sigmaSqz_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
  }else{
    samples_sigmaSqz_r = PROTECT(Rf_allocMatrix(REALSXP, r, nSamples)); nProtect++;
  }
  SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, nr, nSamples)); nProtect++;
  double *Beta = REAL(samples_beta_r);
  double *SigmaSq = REAL(samples_sigmaSq_r);
  double *SigmaSqz = REAL(samples_sigmaSqz_r);
  double *Z = REAL(samples_z_r);

  // Draws are processed in blocks of nb. For each draw s (in order): sigmaSq_s, then beta_s, then the nr standard
  // normal variates of u_s and the n of d_s, so the random-number stream does not depend on the block size.
  const int nbMax = 256;
  int nb = 0, rnb = 0;
  double *Res = (double *) R_alloc((size_t) n * nbMax, sizeof(double));                    // y - X*beta - XTildeBlk*u - d, then w
  double *Tm = (double *) R_alloc((size_t) nr * nbMax, sizeof(double));                    // Rcal*t(XTildeBlk)*w (layout of Z)
  double *sd = (double *) R_alloc(nbMax, sizeof(double));
  double *sdz = (double *) R_alloc(r, sizeof(double));                                     // 1/sqrt(deltasq_q)
  for(q = 0; q < r; q++){
    sdz[q] = 1.0 / sqrt(deltasq[q]);
  }
  double *Zb = NULL, *Tb = NULL, *res_b = NULL, *z_b = NULL, *beta_s = NULL;

  GetRNGstate();

  for(int s0 = 0; s0 < nSamples; s0 += nbMax){

    nb = std::min(nbMax, nSamples - s0);
    rnb = r * nb;
    Zb = Z + (size_t) s0 * nr;

    for(b = 0; b < nb; b++){
      s = s0 + b;
      dtemp = 1.0 / rgamma(sigmaSqIGaPost, 1.0 / sigmaSqIGbPost);
      SigmaSq[s] = dtemp;
      if(nR == 1){
        SigmaSqz[s] = dtemp / deltasq[0];
      }else{
        for(q = 0; q < r; q++){
          SigmaSqz[(size_t) s * r + q] = dtemp / deltasq[q];
        }
      }
      sd[b] = sqrt(dtemp);
      beta_s = Beta + (size_t) s * p;
      for(j = 0; j < p; j++){
        beta_s[j] = rnorm(tmp_p1[j], sd[b]);
      }
      F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, beta_s, &incOne FCONE FCONE FCONE);     // beta_s ~ N(betahat, sigmaSq_s*B)
      z_b = Zb + (size_t) b * nr;
      for(k = 0; k < nr; k++){
        z_b[k] = norm_rand();
      }
      res_b = Res + (size_t) b * n;
      for(i = 0; i < n; i++){
        res_b[i] = Y[i] - sd[b] * norm_rand();                                                       // Y - d, d ~ N(0, sigmaSq_s*I)
      }
    }

    // u ~ N(0, sigmaSq*Rcal): u_q = (sqrt(sigmaSq)/sqrt(deltasq_q))*chol(R_q)*e_q
    if(nR == 1){
      F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &n, &rnb, &one, Rbuf, &n, Zb, &n FCONE FCONE FCONE FCONE);
    }else{
      for(q = 0; q < r; q++){
        F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &n, &nb, &one, Rbuf + nn * q, &n, Zb + (size_t) q * n, &nr FCONE FCONE FCONE FCONE);
      }
    }
    for(b = 0; b < nb; b++){
      for(q = 0; q < r; q++){
        dtemp = sd[b] * sdz[q];
        F77_NAME(dscal)(&n, &dtemp, Zb + (size_t) b * nr + (size_t) q * n, &incOne);
      }
    }

    // Res = Y - d - X*beta - XTildeBlk*u
    F77_NAME(dgemm)(ntran, ntran, &n, &nb, &p, &negOne, X, &n, Beta + (size_t) s0 * p, &p, &one, Res, &n FCONE FCONE);
    for(b = 0; b < nb; b++){
      res_b = Res + (size_t) b * n;
      z_b = Zb + (size_t) b * nr;
      for(q = 0; q < r; q++){
        xq = XTilde + (size_t) n * q;
        const double *u_q = z_b + (size_t) q * n;
        for(i = 0; i < n; i++){
          res_b[i] -= xq[i] * u_q[i];
        }
      }
    }

    // w = inv(Vy)*Res
    solveUt(F, nb, Res, n);
    solveU(F, nb, Res, n);

    // Tm_q = R_q*D_q*w = chol(R_q)*t(chol(R_q))*D_q*w; Tm has the layout of Z (nr x nb), so if R is shared all r
    // blocks are multiplied at once as an n x (r*nb) matrix
    for(b = 0; b < nb; b++){
      res_b = Res + (size_t) b * n;
      Tb = Tm + (size_t) b * nr;
      for(q = 0; q < r; q++){
        xq = XTilde + (size_t) n * q;
        for(i = 0; i < n; i++){
          Tb[(size_t) q * n + i] = xq[i] * res_b[i];
        }
      }
    }
    if(nR == 1){
      F77_NAME(dtrmm)(lside, lower, ytran, nUnit, &n, &rnb, &one, Rbuf, &n, Tm, &n FCONE FCONE FCONE FCONE);
      F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &n, &rnb, &one, Rbuf, &n, Tm, &n FCONE FCONE FCONE FCONE);
    }else{
      for(q = 0; q < r; q++){
        F77_NAME(dtrmm)(lside, lower, ytran, nUnit, &n, &nb, &one, Rbuf + nn * q, &n, Tm + (size_t) q * n, &nr FCONE FCONE FCONE FCONE);
        F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &n, &nb, &one, Rbuf + nn * q, &n, Tm + (size_t) q * n, &nr FCONE FCONE FCONE FCONE);
      }
    }

    // z = u + Rcal*t(XTildeBlk)*w
    for(b = 0; b < nb; b++){
      z_b = Zb + (size_t) b * nr;
      Tb = Tm + (size_t) b * nr;
      for(q = 0; q < r; q++){
        dtemp = 1.0 / deltasq[q];
        F77_NAME(daxpy)(&n, &dtemp, Tb + (size_t) q * n, &incOne, z_b + (size_t) q * n, &incOne);
      }
    }

  }

  PutRNGstate();

  /*****************************************
   PSIS leave-one-out predictive densities
   *****************************************/
  if(loopd_psis){

    int psis_L = psis_tail_length(nSamples);
    int psis_M = psis_gpd_grid_length(psis_L);
    double *ll_i = (double *) R_alloc(nSamples, sizeof(double));
    double *lw_i = (double *) R_alloc(nSamples, sizeof(double));
    int *idx_i = (int *) R_alloc(nSamples, sizeof(int));
    double *xtail_i = (double *) R_alloc(psis_L, sizeof(double));
    double *theta_gpd = (double *) R_alloc(psis_M, sizeof(double));
    double *ltheta_gpd = (double *) R_alloc(psis_M, sizeof(double));
    double theta_i = 0.0;
    double *z_s = NULL;

    loopd_k_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

    for(i = 0; i < n; i++){
      for(s = 0; s < nSamples; s++){
        beta_s = Beta + (size_t) s * p;
        z_s = Z + (size_t) s * nr;
        theta_i = 0.0;
        for(j = 0; j < p; j++){
          theta_i += X[(size_t) j * n + i] * beta_s[j];
        }
        for(q = 0; q < r; q++){
          theta_i += XTilde[(size_t) q * n + i] * z_s[(size_t) q * n + i];
        }
        ll_i[s] = Rf_dnorm4(Y[i], theta_i, sqrt(SigmaSq[s]), 1);
      }
      psis_loo(ll_i, nSamples, psis_L, lw_i, idx_i, xtail_i, theta_gpd, ltheta_gpd,
               &REAL(loopd_out_r)[i], &REAL(loopd_k_r)[i]);
    }

  }

  /*****************************************
   Return object
   *****************************************/
  int nOut = 4 + (loopd ? 1 : 0) + (loopd_psis ? 1 : 0);
  SEXP result_r = PROTECT(Rf_allocVector(VECSXP, nOut)); nProtect++;
  SEXP resultName_r = PROTECT(Rf_allocVector(STRSXP, nOut)); nProtect++;
  SET_VECTOR_ELT(result_r, 0, samples_beta_r);     SET_STRING_ELT(resultName_r, 0, Rf_mkChar("beta"));
  SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);  SET_STRING_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));
  SET_VECTOR_ELT(result_r, 2, samples_sigmaSqz_r); SET_STRING_ELT(resultName_r, 2, Rf_mkChar("sigmaSq.z"));
  SET_VECTOR_ELT(result_r, 3, samples_z_r);        SET_STRING_ELT(resultName_r, 3, Rf_mkChar("z"));
  if(loopd){
    SET_VECTOR_ELT(result_r, 4, loopd_out_r);      SET_STRING_ELT(resultName_r, 4, Rf_mkChar("loopd"));
  }
  if(loopd_psis){
    SET_VECTOR_ELT(result_r, 5, loopd_k_r);        SET_STRING_ELT(resultName_r, 5, Rf_mkChar("loopd.pareto_k"));
  }
  Rf_setAttrib(result_r, R_NamesSymbol, resultName_r);

  // diagnostics, one row per distinct correlation matrix: smallest relative Cholesky pivot of chol(R_k) and chol(Vy),
  // and the off-diagonal range of R_k
  double *diagPivot = (double *) R_alloc(nR, sizeof(double));
  for(k = 0; k < nR; k++){
    diagPivot[k] = fmin2(pivR[k], pivVy);
  }
  result_r = PROTECT(appendDiagnosticsRows(result_r, diagPivot, minCor, maxCor, nR)); nProtect++;

  UNPROTECT(nProtect);
  return result_r;

}

// Reads the common arguments, builds and factorizes the correlation matrices once, and fits the G candidate models
// given by the columns of deltasq (nD x G, nD = r for 'independent' and 1 for 'independent.shared').
static SEXP stvcLM_run(SEXP Y_r, SEXP X_r, SEXP XTilde_r, SEXP n_r, SEXP p_r, SEXP r_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r, SEXP processType_r,
                       SEXP phi_s_r, SEXP phi_t_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                       SEXP deltasq_r, SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r, int verbose){

  int g, q;
  const int incOne = 1;

  double *Y = REAL(Y_r);
  double *X = REAL(X_r);
  double *XTilde = REAL(XTilde_r);
  int n = INTEGER(n_r)[0];
  int p = INTEGER(p_r)[0];
  int r = INTEGER(r_r)[0];
  int pp = p * p;
  double *coords_sp = REAL(sp_coords_r);
  double *coords_tm = REAL(time_coords_r);
  std::string corfn = CHAR(STRING_ELT(corfn_r, 0));
  if(corfn != "gneiting-decay"){
    Rf_error("c++ error: cor.fn must be 'gneiting-decay'.");
  }
  std::string processType = CHAR(STRING_ELT(processType_r, 0));
  int nR = 0;
  if(processType == "independent"){
    nR = r;
  }else if(processType == "independent.shared"){
    nR = 1;
  }else{
    Rf_error("c++ error: process.type must be 'independent' or 'independent.shared'.");
  }
  if(Rf_length(phi_s_r) != nR || Rf_length(phi_t_r) != nR){
    Rf_error("c++ error: phi_s and phi_t must be of length %i.", nR);
  }
  double *phi_s = REAL(phi_s_r);
  double *phi_t = REAL(phi_t_r);

  std::string betaPrior = CHAR(STRING_ELT(betaPrior_r, 0));
  double *betaMu = NULL, *betaV = NULL;
  if(betaPrior == "normal"){
    betaMu = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
    betaV = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaV, &incOne);
  }
  double sigmaSqIGa = REAL(sigmaSqIG_r)[0];
  double sigmaSqIGb = REAL(sigmaSqIG_r)[1];

  int nD = nR;
  if(Rf_length(deltasq_r) % nD != 0 || Rf_length(deltasq_r) == 0){
    Rf_error("c++ error: the length of deltasq must be a multiple of %i.", nD);
  }
  int G = Rf_length(deltasq_r) / nD;
  double *deltasqAll = REAL(deltasq_r);

  int nSamples = INTEGER(nSamples_r)[0];
  int loopd = INTEGER(loopd_r)[0];
  std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));

  if(verbose){
    Rprintf("----------------------------------------\n");
    Rprintf("\tModel description\n");
    Rprintf("----------------------------------------\n");
    Rprintf("Model fit with %i observations.\n\n", n);
    Rprintf("Number of covariates %i (including intercept).\n", p);
    Rprintf("Number of covariates with spatial-temporally varying coefficients %i.\n\n", r);
    Rprintf("Using the %s spatial-temporal correlation function.\n", corfn.c_str());
    Rprintf("Process type: %s.\n\n", processType.c_str());
    Rprintf("Priors:\n");
    if(betaPrior == "flat"){
      Rprintf("\tbeta flat.\n");
    }else{
      Rprintf("\tbeta: Gaussian\n");
      Rprintf("\tmu:"); printVec(betaMu, p);
      Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
      Rprintf("\n");
    }
    if(sigmaSqIGb == 0.0){
      Rprintf("\tsigma.sq: flat, proportional to 1/sigma.sq.\n\n");
    }else{
      Rprintf("\tsigma.sq: Inverse-Gamma\n\tshape = %.2f, scale = %.2f.\n\n", sigmaSqIGa, sigmaSqIGb);
    }
    Rprintf("Spatial-temporal process parameters:\n");
    for(q = 0; q < nR; q++){
      if(nR > 1){
        Rprintf("\tprocess %i: phi_s = %.2f, phi_t = %.2f, noise-to-spatial variance ratio = %.2f.\n",
                q + 1, phi_s[q], phi_t[q], deltasqAll[q]);
      }else{
        Rprintf("\tphi_s = %.2f, phi_t = %.2f, noise-to-spatial variance ratio = %.2f.\n",
                phi_s[q], phi_t[q], deltasqAll[q]);
      }
    }
    Rprintf("\nNumber of posterior samples = %i.\n\n", nSamples);
    if(loopd){
      Rprintf("LOO-PD calculation method = %s.\n", loopd_method.c_str());
    }
    Rprintf("----------------------------------------\n");
  }

  // correlation matrices and their Cholesky factors, shared by the G candidate models
  double *Rbuf = (double *) R_alloc((size_t) n * n * nR, sizeof(double));
  double *pivR = (double *) R_alloc(nR, sizeof(double));
  double *minCor = (double *) R_alloc(nR, sizeof(double));
  double *maxCor = (double *) R_alloc(nR, sizeof(double));
  stvcLM_corChol(n, nR, coords_sp, coords_tm, phi_s, phi_t, corfn, Rbuf, pivR, minCor, maxCor);

  double *deltasq = (double *) R_alloc(r, sizeof(double));
  SEXP result_r = PROTECT(Rf_allocVector(VECSXP, G));

  for(g = 0; g < G; g++){

    const void *vmax = vmaxget();                                                  // release the per-model R_alloc memory below
    for(q = 0; q < r; q++){
      deltasq[q] = deltasqAll[(size_t) g * nD + (nD == 1 ? 0 : q)];
    }
    SET_VECTOR_ELT(result_r, g, stvcLMexact_fit(Y, X, XTilde, n, p, r, nR, Rbuf, pivR, minCor, maxCor,
                                                betaPrior, betaMu, betaV, sigmaSqIGa, sigmaSqIGb, deltasq,
                                                nSamples, loopd, loopd_method));
    vmaxset(vmax);

    R_CheckUserInterrupt();

  }

  UNPROTECT(1);
  return result_r;

}

extern "C" {

  // Single model: deltasq has length r ('independent') or 1 ('independent.shared'); returns the fit.
  SEXP stvcLMexact(SEXP Y_r, SEXP X_r, SEXP XTilde_r, SEXP n_r, SEXP p_r, SEXP r_r,
                   SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r, SEXP processType_r,
                   SEXP phi_s_r, SEXP phi_t_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                   SEXP deltasq_r, SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r, SEXP verbose_r){

    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    int nD = (processType == "independent") ? INTEGER(r_r)[0] : 1;
    if(Rf_length(deltasq_r) != nD){
      Rf_error("c++ error: deltasq must be of length %i.", nD);
    }
    SEXP out_r = PROTECT(stvcLM_run(Y_r, X_r, XTilde_r, n_r, p_r, r_r, sp_coords_r, time_coords_r, corfn_r,
                                    processType_r, phi_s_r, phi_t_r, betaPrior_r, betaNorm_r, sigmaSqIG_r,
                                    deltasq_r, nSamples_r, loopd_r, loopd_method_r, INTEGER(verbose_r)[0]));
    UNPROTECT(1);
    return VECTOR_ELT(out_r, 0);

  }

  // Candidate models sharing (phi_s, phi_t): deltasq is an nD x G matrix (one column per candidate model); the
  // correlation matrices and their Cholesky factors are computed once. Returns a list of G fits.
  SEXP stvcLMexactGrid(SEXP Y_r, SEXP X_r, SEXP XTilde_r, SEXP n_r, SEXP p_r, SEXP r_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r, SEXP processType_r,
                       SEXP phi_s_r, SEXP phi_t_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                       SEXP deltasq_r, SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r){

    return stvcLM_run(Y_r, X_r, XTilde_r, n_r, p_r, r_r, sp_coords_r, time_coords_r, corfn_r,
                      processType_r, phi_s_r, phi_t_r, betaPrior_r, betaNorm_r, sigmaSqIG_r,
                      deltasq_r, nSamples_r, loopd_r, loopd_method_r, 0);

  }

}
