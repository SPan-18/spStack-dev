#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
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

extern "C" {

    SEXP predict_spGLM(SEXP n_r, SEXP n_pred_r, SEXP p_r, SEXP family_r, SEXP nBinom_new_r,
                       SEXP X_new_r, SEXP sp_coords_r, SEXP sp_coords_new_r,
                       SEXP corfn_r, SEXP phi_r, SEXP nu_r, SEXP nSamples_r,
                       SEXP beta_samps_r, SEXP z_samps_r, SEXP sigmaSq_z_samps_r, SEXP joint_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, s, info, nProtect = 0;
    char const *lower = "L";
    char const *nunit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int n_pred = INTEGER(n_pred_r)[0];
    int n_predn_pred = n_pred * n_pred;
    int nn_pred = n * n_pred;
    int p = INTEGER(p_r)[0];
    int joint = INTEGER(joint_r)[0];
    double *X_new = REAL(X_new_r);
    int *nBinom_new = INTEGER(nBinom_new_r);

    double *zSamps = REAL(z_samps_r);
    double *betaSamps = REAL(beta_samps_r);
    double *sigmaSqzSamps = REAL(sigmaSq_z_samps_r);

    double *coords_sp = REAL(sp_coords_r);
    double *coords_sp_new = REAL(sp_coords_new_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    std::string family = CHAR(STRING_ELT(family_r, 0));
    const char *family_poisson = "poisson";
    const char *family_binary = "binary";
    const char *family_binomial = "binomial";

    // spatial process parameters
    double phi = REAL(phi_r)[0];
    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];
    double thetasp[2] = {phi, nu};

    // Given a posterior draw (beta, z, sigmaSqz):
    //   z.pred | z ~ N(t(C)*inv(R)*z, sigmaSqz*(R_new - t(C)*inv(R)*C)),  mu.pred = g^{-1}(X_new*beta + z.pred),
    //   y.pred ~ Poisson(mu.pred) or Binomial(nBinom_new, mu.pred),
    // where R = cor(observed), C = cor(observed, new) and R_new = cor(new) (joint) or its diagonal of ones (pointwise).

    // chol(R), built in place
    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);
    spCorFull2(n, 2, coords_sp, thetasp, corfn, cholVz);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE);
    if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial correlation matrix failed (info = %i).\n", info);}

    // Cz = cholinv(R)*C
    double *Cz = (double *) R_alloc(nn_pred, sizeof(double)); zeros(Cz, nn_pred);
    spCorCross(n, n_pred, 2, coords_sp, coords_sp_new, thetasp, corfn, Cz);
    F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);

    // conditional covariance of z.pred given z (per unit sigmaSqz)
    double *z_pred_cov = NULL;
    if(joint){
      // z_pred_cov = chol(R_new - t(Cz)*Cz), formed in place (lower triangle)
      z_pred_cov = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(z_pred_cov, n_predn_pred);
      spCorFull2(n_pred, 2, coords_sp_new, thetasp, corfn, z_pred_cov);                                                    // z_pred_cov = R_new
      F77_NAME(dsyrk)(lower, ytran, &n_pred, &n, &negOne, Cz, &n, &one, z_pred_cov, &n_pred FCONE FCONE);                 // z_pred_cov = R_new - t(Cz)*Cz
      F77_NAME(dpotrf)(lower, &n_pred, z_pred_cov, &n_pred, &info FCONE);
      if(info != 0){Rf_error("c++ error: Cholesky factorization of the conditional covariance of z.pred failed (info = %i); check for prediction locations that coincide with each other or with observed locations.\n", info);}
      mkLT(z_pred_cov, n_pred);
    }else{
      // pointwise conditional variances 1 - ||Cz[, i]||^2, clamped at 0 against rounding
      z_pred_cov = (double *) R_alloc(n_pred, sizeof(double)); zeros(z_pred_cov, n_pred);
      for(i = 0; i < n_pred; i++){
        z_pred_cov[i] = fmax2(1.0 - F77_CALL(ddot)(&n, &Cz[i * n], &incOne, &Cz[i * n], &incOne), 0.0);
      }
    }

    // posterior predictive samples of z and y
    SEXP samples_predz_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predmu_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predy_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    double *predz = REAL(samples_predz_r);
    double *predmu = REAL(samples_predmu_r);
    double *predy = REAL(samples_predy_r);

    // Draws are processed in blocks of nBlock: the conditional means t(Cz)*cholinv(R)*Z of the whole block are
    // found with level-3 BLAS; the random variates are then drawn draw by draw, in the same order as before
    // (for each draw, the n_pred variates of z.pred and then the n_pred of y.pred, whose number depends on mu.pred).
    const int nBlockMax = 64;
    int nBlock = 0, nzBlock = 0, b = 0;
    R_xlen_t offset = 0;
    double dtemp1 = 0.0;
    double *zBlock = (double *) R_alloc((size_t) n * nBlockMax, sizeof(double));                    // n x nBlock work matrix
    double *noise = (double *) R_alloc(n_pred, sizeof(double)); zeros(noise, n_pred);
    double *zp = NULL;

    GetRNGstate();

    for(s = 0; s < nSamples; s += nBlockMax){

      nBlock = std::min(nBlockMax, nSamples - s);
      nzBlock = n * nBlock;
      offset = (R_xlen_t) s * n_pred;

      // conditional means of the block: predz[, s:(s+nBlock)] = t(Cz)*cholinv(R)*Z[, s:(s+nBlock)]
      F77_NAME(dcopy)(&nzBlock, &zSamps[(R_xlen_t) s * n], &incOne, zBlock, &incOne);
      F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &nBlock, &one, cholVz, &n, zBlock, &n FCONE FCONE FCONE FCONE);
      F77_NAME(dgemm)(ytran, ntran, &n_pred, &nBlock, &n, &one, Cz, &n, zBlock, &n, &zero, &predz[offset], &n_pred FCONE FCONE);

      for(b = 0; b < nBlock; b++){

        zp = &predz[offset + (R_xlen_t) b * n_pred];

        // z.pred = mean + noise
        if(joint){
          for(i = 0; i < n_pred; i++){
            noise[i] = rnorm(0.0, sqrt(sigmaSqzSamps[s + b]));
          }
          F77_NAME(dtrmv)(lower, ntran, nunit, &n_pred, z_pred_cov, &n_pred, noise, &incOne FCONE FCONE FCONE);   // noise = chol(z_pred_cov)*noise
          F77_NAME(daxpy)(&n_pred, &one, noise, &incOne, zp, &incOne);
        }else{
          for(i = 0; i < n_pred; i++){
            zp[i] += rnorm(0.0, sqrt(sigmaSqzSamps[s + b] * z_pred_cov[i]));
          }
        }

        // natural parameter X_new*beta + z.pred, then mu.pred and y.pred
        F77_NAME(dcopy)(&n_pred, zp, &incOne, noise, &incOne);
        F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, &betaSamps[(R_xlen_t) (s + b) * p], &incOne, &one, noise, &incOne FCONE);
        for(i = 0; i < n_pred; i++){
          if(family == family_poisson){
            dtemp1 = exp(noise[i]);
            predmu[offset + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offset + (R_xlen_t) b * n_pred + i] = rpois(dtemp1);
          }else if(family == family_binary){
            dtemp1 = inverse_logit(noise[i]);
            predmu[offset + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offset + (R_xlen_t) b * n_pred + i] = rbinom(1, dtemp1);
          }else if(family == family_binomial){
            dtemp1 = inverse_logit(noise[i]);
            predmu[offset + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offset + (R_xlen_t) b * n_pred + i] = rbinom(nBinom_new[i], dtemp1);
          }
        }

      }

    }

    PutRNGstate();

    // make return object
    SEXP result_r, resultName_r;

    // make return object for posterior samples
    int nResultListObjs = 3;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // posterior predictive samples of spatial-temporal process z
    SET_VECTOR_ELT(result_r, 0, samples_predz_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("z.pred"));

    // posterior predictive samples of the canonical parameter mu
    SET_VECTOR_ELT(result_r, 1, samples_predmu_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("mu.pred"));

    // posterior predictive samples of the response variable y
    SET_VECTOR_ELT(result_r, 2, samples_predy_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("y.pred"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

    }

}