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

    SEXP predict_spLM(SEXP n_r, SEXP n_pred_r, SEXP p_r,
                      SEXP X_new_r, SEXP sp_coords_r, SEXP sp_coords_new_r,
                      SEXP corfn_r, SEXP phi_r, SEXP nu_r, SEXP deltasq_r,
                      SEXP beta_samps_r, SEXP z_samps_r, SEXP sigmaSq_z_samps_r,
                      SEXP nSamples_r, SEXP joint_r){

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

    double *zSamps = REAL(z_samps_r);
    double *betaSamps = REAL(beta_samps_r);
    double *sigmaSqzSamps = REAL(sigmaSq_z_samps_r);

    double *coords_sp = REAL(sp_coords_r);
    double *coords_sp_new = REAL(sp_coords_new_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // spatial process parameters
    double phi = REAL(phi_r)[0];
    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }
    double deltasq = REAL(deltasq_r)[0];

    // Create correlation and cross-correlation matrices
    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
    double *Cz = (double *) R_alloc(nn_pred, sizeof(double)); zeros(Cz, nn_pred);
    double *Vz_new = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(Vz_new, n_predn_pred);
    double *thetasp = (double *) R_alloc(2, sizeof(double));

    //construct covariance matrix (full)
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull2(n, 2, coords_sp, thetasp, corfn, Vz);
    spCorCross(n, n_pred, 2, coords_sp, coords_sp_new, thetasp, corfn, Cz);
    spCorFull2(n_pred, 2, coords_sp_new, thetasp, corfn, Vz_new);

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];

    // posterior predictive samples of z and y
    SEXP samples_predz_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predmu_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predy_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;

    // Set up pre-processing matrices etc.
    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);  // chol(Vz)
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
    mkLT(cholVz, n);

    // Cz = cholinv(Vz)*Cz
    F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);

    double *z_pred_cov = NULL;  // define NULL pointer for z_pred_cov
    double *z_pred_mu = (double *) R_alloc(n_pred, sizeof(double)); zeros(z_pred_mu, n_pred);  // n_predx1 vector z_pred_mu
    double *beta_s = (double *) R_alloc(p, sizeof(double)); zeros(beta_s, p);
    double *z_s = (double *) R_alloc(n, sizeof(double)); zeros(z_s, n);
    double *z_pred_s = (double *) R_alloc(n_pred, sizeof(double)); zeros(z_pred_s, n_pred);      // n_predx1 vector z_pred_s
    double *tmp_n_pred = (double *) R_alloc(n_pred, sizeof(double)); zeros(tmp_n_pred, n_pred);  // n_predx1 vector tmp_n_pred
    double dtemp1 = 0.0;

    if(joint){

        z_pred_cov = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(z_pred_cov, n_predn_pred);  // n_predxn_pred matrix z_pred_cov

        // first find the Schur complement RTilde - t(C)*inv(R)*C
        // goal: z_pred_cov = chol(t(Cz)*inv(Vz)*Cz)
        F77_NAME(dgemm)(ytran, ntran, &n_pred, &n_pred, &n, &one, Cz, &n, Cz, &n, &zero, z_pred_cov, &n_pred FCONE FCONE);
        F77_NAME(daxpy)(&n_predn_pred, &negOne, Vz_new, &incOne, z_pred_cov, &incOne);

        // z_pred_cov = Vz_new - t(Cz)*inv(Vz)*Cz
        F77_NAME(dscal)(&n_predn_pred, &negOne, z_pred_cov, &incOne);
        F77_NAME(dpotrf)(lower, &n_pred, z_pred_cov, &n_pred, &info FCONE); if(info != 0){perror("c++ error: z_pred_cov dpotrf failed\n");}
        mkLT(z_pred_cov, n_pred);

        for(s = 0; s < nSamples; s++){

            // copy posterior samples of beta and z
            F77_NAME(dcopy)(&p, &betaSamps[s * p], &incOne, beta_s, &incOne);
            F77_NAME(dcopy)(&n, &zSamps[s * n], &incOne, z_s, &incOne);

            // find posterior predictive mean z_pred_mu = t(Cz)*inv(Vz)*z_s
            F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVz, &n, z_s, &incOne FCONE FCONE FCONE);
            F77_NAME(dgemv)(ytran, &n, &n_pred, &one, Cz, &n, z_s, &incOne, &zero, z_pred_mu, &incOne FCONE);

            for(i = 0; i < n_pred; i++){
                tmp_n_pred[i] = rnorm(0.0, sqrt(sigmaSqzSamps[s]));
            }
            F77_NAME(dgemv)(ntran, &n_pred, &n_pred, &one, z_pred_cov, &n_pred, tmp_n_pred, &incOne, &zero, z_pred_s, &incOne FCONE);
            F77_NAME(daxpy)(&n_pred, &one, z_pred_mu, &incOne, z_pred_s, &incOne);
            F77_NAME(dcopy)(&n_pred, &z_pred_s[0], &incOne, &REAL(samples_predz_r)[s*n_pred], &incOne);

            // Find natural parameter X*beta + z; dgemv: z_pred_s = z_pred_s + X_new * beta_s
            F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, beta_s, &incOne, &one, z_pred_s, &incOne FCONE);

            // Sample from predictive distribution
            for(i = 0; i < n_pred; i++){

                dtemp1 = z_pred_s[i];
                REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                REAL(samples_predy_r)[s * n_pred + i] = rnorm(dtemp1, sqrt(deltasq * sigmaSqzSamps[s]));

            }

        }
        // End of joint prediction set-up

    }else{

        z_pred_cov = (double *) R_alloc(n_pred, sizeof(double)); zeros(z_pred_cov, n_pred);  // n_predx1 vector z_pred_cov
        for(i = 0; i < n_pred; i++){
            z_pred_cov[i] = 1.0 - F77_CALL(ddot)(&n, &Cz[i * n], &incOne, &Cz[i * n], &incOne);
        }

        for(s = 0; s < nSamples; s++){

            // copy posterior samples of beta and z
            F77_NAME(dcopy)(&p, &betaSamps[s * p], &incOne, beta_s, &incOne);
            F77_NAME(dcopy)(&n, &zSamps[s * n], &incOne, z_s, &incOne);

            // find posterior predictive mean z_pred_mu = t(Cz)*inv(Vz)*z_s
            F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVz, &n, z_s, &incOne FCONE FCONE FCONE);
            F77_NAME(dgemv)(ytran, &n, &n_pred, &one, Cz, &n, z_s, &incOne, &zero, z_pred_mu, &incOne FCONE);

            for(i = 0; i < n_pred; i++){
                z_pred_s[i] = rnorm(0.0, sqrt(sigmaSqzSamps[s] * z_pred_cov[i]));
            }
            F77_NAME(daxpy)(&n_pred, &one, z_pred_mu, &incOne, z_pred_s, &incOne);
            F77_NAME(dcopy)(&n_pred, &z_pred_s[0], &incOne, &REAL(samples_predz_r)[s*n_pred], &incOne);

            // Find natural parameter X*beta + z; dgemv: z_pred_s = z_pred_s + X_new * beta_s
            F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, beta_s, &incOne, &one, z_pred_s, &incOne FCONE);

            // Sample from predictive distribution
            // Sample from predictive distribution
            for(i = 0; i < n_pred; i++){

                dtemp1 = z_pred_s[i];
                REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                REAL(samples_predy_r)[s * n_pred + i] = rnorm(dtemp1, sqrt(deltasq * sigmaSqzSamps[s]));

            }

        }
        // End of point-wise prediction set-up

    }

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