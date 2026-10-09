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

  SEXP predict_stvcGLM(SEXP n_r, SEXP n_pred_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_new_r,
                       SEXP X_new_r, SEXP XTilde_new_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP sp_coords_new_r, SEXP time_coords_new_r,
                       SEXP processType_r, SEXP corfn_r, SEXP phi_s_r, SEXP phi_t_r, SEXP nSamples_r,
                       SEXP beta_samps_r, SEXP z_samps_r, SEXP z_scale_samps_r, SEXP joint_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, k, s, info, nProtect = 0;
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
    int r = INTEGER(r_r)[0];
    int rr = r * r;
    int nr = n * r;
    int n_predr = n_pred * r;
    int joint = INTEGER(joint_r)[0];
    double *X_new = REAL(X_new_r);
    double *XTilde_new = REAL(XTilde_new_r);
    int *nBinom_new = INTEGER(nBinom_new_r);

    double *zSamps = REAL(z_samps_r);
    double *betaSamps = REAL(beta_samps_r);
    double *zScaleSamps = REAL(z_scale_samps_r);

    double *coords_sp = REAL(sp_coords_r);
    double *coords_sp_new = REAL(sp_coords_new_r);
    double *coords_tm = REAL(time_coords_r);
    double *coords_tm_new = REAL(time_coords_new_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    std::string family = CHAR(STRING_ELT(family_r, 0));
    const char *family_poisson = "poisson";
    const char *family_binary = "binary";
    const char *family_binomial = "binomial";

    // supported spatial-temporal process models
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    if(processType != "independent.shared" && processType != "independent" && processType != "multivariate"){
      Rf_error("c++ error: process.type must be one of 'independent', 'independent.shared' or 'multivariate'.");
    }
    if(corfn != "gneiting-decay"){
      Rf_error("c++ error: cor.fn must be 'gneiting-decay'.");
    }
    int nCov = (processType == "independent") ? r : 1;                   // number of distinct correlation functions

    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double thetaspt[2] = {0.0, 0.0};
    F77_NAME(dcopy)(&nCov, REAL(phi_s_r), &incOne, phi_s_vec, &incOne);
    F77_NAME(dcopy)(&nCov, REAL(phi_t_r), &incOne, phi_t_vec, &incOne);

    // Given a posterior draw (beta, z, scale), for each correlation function k (one shared, or r independent):
    //   z.pred_k | z_k ~ N(t(C_k)*inv(R_k)*z_k, scale*(R_new_k - t(C_k)*inv(R_k)*C_k)) (multivariate: matrix normal
    //   with column covariance Sigma), natural parameter X_new*beta + XTilde_new*z.pred, and y.pred from the family.
    // For each k: chol(R_k) built in place, Cz_k = cholinv(R_k)*C_k, and chol(R_new_k - t(Cz_k)*Cz_k) (joint) or the
    // diagonal 1 - ||Cz_k[, i]||^2 clamped at 0 (pointwise).
    double *cholVz = (double *) R_alloc((size_t) nn * nCov, sizeof(double)); zeros(cholVz, nn * nCov);
    double *Cz = (double *) R_alloc((size_t) nn_pred * nCov, sizeof(double)); zeros(Cz, nn_pred * nCov);
    double *z_pred_cov = NULL;
    if(joint){
      z_pred_cov = (double *) R_alloc((size_t) n_predn_pred * nCov, sizeof(double)); zeros(z_pred_cov, n_predn_pred * nCov);
    }else{
      z_pred_cov = (double *) R_alloc((size_t) n_pred * nCov, sizeof(double)); zeros(z_pred_cov, n_pred * nCov);
    }

    for(k = 0; k < nCov; k++){
      thetaspt[0] = phi_s_vec[k];
      thetaspt[1] = phi_t_vec[k];
      sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, &cholVz[nn * k]);
      F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE);
      if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i).\n", info);}
      mkLT(&cholVz[nn * k], n);
      sptCorCross(n, n_pred, 2, coords_sp, coords_tm, coords_sp_new, coords_tm_new, thetaspt, corfn, &Cz[nn_pred * k]);
      F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, &cholVz[nn * k], &n, &Cz[nn_pred * k], &n FCONE FCONE FCONE FCONE);  // Cz = cholinv(R)*C
      if(joint){
        sptCorFull(n_pred, 2, coords_sp_new, coords_tm_new, thetaspt, corfn, &z_pred_cov[n_predn_pred * k]);                        // R_new
        F77_NAME(dsyrk)(lower, ytran, &n_pred, &n, &negOne, &Cz[nn_pred * k], &n, &one, &z_pred_cov[n_predn_pred * k], &n_pred FCONE FCONE);
        F77_NAME(dpotrf)(lower, &n_pred, &z_pred_cov[n_predn_pred * k], &n_pred, &info FCONE);
        if(info != 0){Rf_error("c++ error: Cholesky factorization of the conditional covariance of z.pred failed (info = %i); check for prediction locations that coincide with each other or with observed locations.\n", info);}
        mkLT(&z_pred_cov[n_predn_pred * k], n_pred);
      }else{
        for(i = 0; i < n_pred; i++){
          z_pred_cov[n_pred * k + i] = fmax2(1.0 - F77_CALL(ddot)(&n, &Cz[nn_pred * k + i * n], &incOne, &Cz[nn_pred * k + i * n], &incOne), 0.0);
        }
      }
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];
    // posterior predictive samples of z and y
    SEXP samples_predz_r = PROTECT(Rf_allocMatrix(REALSXP, n_predr, nSamples)); nProtect++;
    SEXP samples_predmu_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predy_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    double *predz = REAL(samples_predz_r);
    double *predmu = REAL(samples_predmu_r);
    double *predy = REAL(samples_predy_r);

    double *zScale_s = (double *) R_alloc(rr, sizeof(double)); zeros(zScale_s, rr);
    double *noise = (double *) R_alloc(n_predr, sizeof(double)); zeros(noise, n_predr);
    double *tmp_n_predr = (double *) R_alloc(n_predr, sizeof(double)); zeros(tmp_n_predr, n_predr);
    double dtemp1 = 0.0;

    // Draws are processed in blocks of nBlock: the conditional means t(Cz_k)*cholinv(R_k)*z_k of the whole block are
    // found with level-3 BLAS; the random variates are then drawn draw by draw, in the same order as before (for each
    // draw, the variates of z.pred and then those of y.pred, whose number depends on mu.pred).
    const int nBlockMax = 64;
    int nBlock = 0, nrBlock = 0, b = 0;
    R_xlen_t offz = 0, offy = 0;
    double *zBlock = (double *) R_alloc((size_t) nr * nBlockMax, sizeof(double));                   // nr x nBlock work matrix
    double *zp = NULL;

    GetRNGstate();

    for(s = 0; s < nSamples; s += nBlockMax){

      nBlock = std::min(nBlockMax, nSamples - s);
      nrBlock = nr * nBlock;
      offz = (R_xlen_t) s * n_predr;
      offy = (R_xlen_t) s * n_pred;

      // conditional means of the block, written into predz
      F77_NAME(dcopy)(&nrBlock, &zSamps[(R_xlen_t) s * nr], &incOne, zBlock, &incOne);
      if(nCov == 1){
        // the r processes share R: treat the block as an n x (r*nBlock) matrix
        int rBlock = r * nBlock;
        F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &rBlock, &one, cholVz, &n, zBlock, &n FCONE FCONE FCONE FCONE);
        F77_NAME(dgemm)(ytran, ntran, &n_pred, &rBlock, &n, &one, Cz, &n, zBlock, &n, &zero, &predz[offz], &n_pred FCONE FCONE);
      }else{
        for(k = 0; k < r; k++){
          F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &nBlock, &one, &cholVz[nn * k], &n, &zBlock[n * k], &nr FCONE FCONE FCONE FCONE);
          F77_NAME(dgemm)(ytran, ntran, &n_pred, &nBlock, &n, &one, &Cz[nn_pred * k], &n, &zBlock[n * k], &nr, &zero, &predz[offz + n_pred * k], &n_predr FCONE FCONE);
        }
      }

      for(b = 0; b < nBlock; b++){

        zp = &predz[offz + (R_xlen_t) b * n_predr];

        // z.pred = mean + noise, with the variates drawn in the original order
        if(processType == "independent.shared"){

          dtemp1 = zScaleSamps[s + b];
          if(joint){
            for(k = 0; k < r; k++){
              for(i = 0; i < n_pred; i++){
                noise[k * n_pred + i] = rnorm(0.0, sqrt(dtemp1));
              }
            }
            F77_NAME(dtrmm)(lside, lower, ntran, nunit, &n_pred, &r, &one, z_pred_cov, &n_pred, noise, &n_pred FCONE FCONE FCONE FCONE);
            F77_NAME(daxpy)(&n_predr, &one, noise, &incOne, zp, &incOne);
          }else{
            for(k = 0; k < r; k++){
              for(i = 0; i < n_pred; i++){
                zp[k * n_pred + i] += rnorm(0.0, sqrt(dtemp1 * z_pred_cov[i]));
              }
            }
          }

        }else if(processType == "independent"){

          for(k = 0; k < r; k++){
            dtemp1 = zScaleSamps[(R_xlen_t) (s + b) * r + k];
            if(joint){
              for(i = 0; i < n_pred; i++){
                noise[i] = rnorm(0.0, sqrt(dtemp1));
              }
              F77_NAME(dtrmv)(lower, ntran, nunit, &n_pred, &z_pred_cov[n_predn_pred * k], &n_pred, noise, &incOne FCONE FCONE FCONE);
              F77_NAME(daxpy)(&n_pred, &one, noise, &incOne, &zp[n_pred * k], &incOne);
            }else{
              for(i = 0; i < n_pred; i++){
                zp[k * n_pred + i] += rnorm(0.0, sqrt(dtemp1 * z_pred_cov[k * n_pred + i]));
              }
            }
          }

        }else if(processType == "multivariate"){

          // column covariance Sigma of the draw
          F77_NAME(dcopy)(&rr, &zScaleSamps[(R_xlen_t) (s + b) * rr], &incOne, zScale_s, &incOne);
          F77_NAME(dpotrf)(lower, &r, zScale_s, &r, &info FCONE);
          if(info != 0){PutRNGstate(); Rf_error("c++ error: Cholesky factorization of a posterior sample of Sigma failed (info = %i).\n", info);}
          mkLT(zScale_s, r);
          for(k = 0; k < r; k++){
            for(i = 0; i < n_pred; i++){
              if(joint){
                noise[k * n_pred + i] = rnorm(0.0, 1.0);
              }else{
                noise[k * n_pred + i] = rnorm(0.0, sqrt(z_pred_cov[i]));
              }
            }
          }
          if(joint){
            F77_NAME(dtrmm)(lside, lower, ntran, nunit, &n_pred, &r, &one, z_pred_cov, &n_pred, noise, &n_pred FCONE FCONE FCONE FCONE);   // chol(z_pred_cov)*noise
          }
          F77_NAME(dgemm)(ntran, ytran, &n_pred, &r, &r, &one, noise, &n_pred, zScale_s, &r, &zero, tmp_n_predr, &n_pred FCONE FCONE);  // *t(chol(Sigma))
          F77_NAME(daxpy)(&n_predr, &one, tmp_n_predr, &incOne, zp, &incOne);

        }

        // natural parameter X_new*beta + XTilde_new*z.pred, then mu.pred and y.pred
        lmulm_XTilde_VC(ntran, n_pred, r, 1, XTilde_new, zp, tmp_n_predr);
        F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, &betaSamps[(R_xlen_t) (s + b) * p], &incOne, &one, tmp_n_predr, &incOne FCONE);
        for(i = 0; i < n_pred; i++){
          if(family == family_poisson){
            dtemp1 = exp(tmp_n_predr[i]);
            predmu[offy + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offy + (R_xlen_t) b * n_pred + i] = rpois(dtemp1);
          }else if(family == family_binary){
            dtemp1 = inverse_logit(tmp_n_predr[i]);
            predmu[offy + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offy + (R_xlen_t) b * n_pred + i] = rbinom(1, dtemp1);
          }else if(family == family_binomial){
            dtemp1 = inverse_logit(tmp_n_predr[i]);
            predmu[offy + (R_xlen_t) b * n_pred + i] = dtemp1;
            predy[offy + (R_xlen_t) b * n_pred + i] = rbinom(nBinom_new[i], dtemp1);
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