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
    int nnr = nn * r;
    int n_predn_predr = n_predn_pred * r;
    int nn_predr = nn_pred * r;
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

    // create spatial-temporal covariance matrices
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    double *Vz = NULL;
    double *Vz_new = NULL;
    double *Cz = NULL;

    if(corfn == "gneiting-decay"){

        if(processType == "independent.shared" || processType == "multivariate"){

            phi_s_vec[0] = REAL(phi_s_r)[0];
            phi_t_vec[0] = REAL(phi_t_r)[0];
            thetaspt[0] = phi_s_vec[0];
            thetaspt[1] = phi_t_vec[0];

            Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
            sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, Vz);

            if(joint){
                Vz_new = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(Vz_new, n_predn_pred);
                sptCorFull(n_pred, 2, coords_sp_new, coords_tm_new, thetaspt, corfn, Vz_new);
            }

            Cz = (double *) R_alloc(nn_pred, sizeof(double)); zeros(Cz, nn_pred);
            sptCorCross(n, n_pred, 2, coords_sp, coords_tm, coords_sp_new, coords_tm_new, thetaspt, corfn, Cz);

        }else if(processType == "independent"){

            F77_NAME(dcopy)(&r, REAL(phi_s_r), &incOne, phi_s_vec, &incOne);
            F77_NAME(dcopy)(&r, REAL(phi_t_r), &incOne, phi_t_vec, &incOne);

            // find r-many correlation/cross-correlation matrices, stacked into a rn^2-dim vector
            Vz = (double *) R_alloc(nnr, sizeof(double)); zeros(Vz, nnr);
            Cz = (double *) R_alloc(nn_predr, sizeof(double)); zeros(Cz, nn_predr);
            for(k = 0; k < r; k++){
                thetaspt[0] = phi_s_vec[k];
                thetaspt[1] = phi_t_vec[k];
                sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, &Vz[nn * k]);
                sptCorCross(n, n_pred, 2, coords_sp, coords_tm, coords_sp_new, coords_tm_new, thetaspt, corfn, &Cz[nn_pred * k]);
            }

            if(joint){
                Vz_new = (double *) R_alloc(n_predn_predr, sizeof(double)); zeros(Vz_new, n_predn_predr);
                for(k = 0; k < r; k++){
                    thetaspt[0] = phi_s_vec[k];
                    thetaspt[1] = phi_t_vec[k];
                    sptCorFull(n_pred, 2, coords_sp_new, coords_tm_new, thetaspt, corfn, &Vz_new[n_predn_pred * k]);
                }
            }

        }
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];
    // posterior predictive samples of z and y
    SEXP samples_predz_r = PROTECT(Rf_allocMatrix(REALSXP, n_predr, nSamples)); nProtect++;
    SEXP samples_predmu_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;
    SEXP samples_predy_r = PROTECT(Rf_allocMatrix(REALSXP, n_pred, nSamples)); nProtect++;

    double *beta_s = (double *) R_alloc(p, sizeof(double)); zeros(beta_s, p);
    double *z_s = (double *) R_alloc(nr, sizeof(double)); zeros(z_s, nr);
    double *zScale_s = NULL;
    if(processType == "independent"){
        zScale_s = (double *) R_alloc(r, sizeof(double)); zeros(zScale_s, r);
    }else if(processType == "multivariate"){
        zScale_s = (double *) R_alloc(rr, sizeof(double)); zeros(zScale_s, rr);
    }
    double dtemp1 = 0.0;

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/

    double *cholVz = NULL;               // define NULL pointer for chol(Vz)
    double *z_pred_mu = (double *) R_alloc(n_predr, sizeof(double)); zeros(z_pred_mu, n_predr);      // n_pred*rx1 vector z_pred_mu
    double *z_pred_s = (double *) R_alloc(n_predr, sizeof(double)); zeros(z_pred_s, n_predr);        // n_pred*rx1 vector z_pred_s
    double *tmp_n_predr = (double *) R_alloc(n_predr, sizeof(double)); zeros(tmp_n_predr, n_predr);  // n_pred*rx1 vector tmp_n_predr

    if(joint){

        // Joint prediction
        double *z_pred_cov = NULL;           // define NULL pointer for z_pred_cov

        // Find Cholesky of Vz, find cholesky of z_pred_cov
        if(processType == "independent.shared" || processType == "multivariate"){

            cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);                              // nxn matrix chol(Vz)
            z_pred_cov = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(z_pred_cov, n_predn_pred);  // n_predxn_pred matrix chol(z_pred_cov)

            F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
            F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
            mkLT(cholVz, n);

            F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);          // Cz = cholinv(Vz)*Cz
            F77_NAME(dgemm)(ytran, ntran, &n_pred, &n_pred, &n, &one, Cz, &n, Cz, &n, &zero, z_pred_cov, &n_pred FCONE FCONE);   // z_pred_cov = t(Cz)*inv(Vz)*Cz
            F77_NAME(daxpy)(&n_predn_pred, &negOne, Vz_new, &incOne, z_pred_cov, &incOne);
            F77_NAME(dscal)(&n_predn_pred, &negOne, z_pred_cov, &incOne);                                                        // z_pred_cov = Vz_new - t(Cz)*inv(Vz)*Cz
            F77_NAME(dpotrf)(lower, &n_pred, z_pred_cov, &n_pred, &info FCONE); if(info != 0){perror("c++ error: z_pred_cov dpotrf failed\n");}
            mkLT(z_pred_cov, n_pred);

        }else if(processType == "independent"){

            cholVz = (double *) R_alloc(nnr, sizeof(double)); zeros(cholVz, nnr);                             // r nxn matrices chol(Vz)
            z_pred_cov = (double *) R_alloc(n_predn_predr, sizeof(double)); zeros(z_pred_cov, n_predn_predr); // r n_predxn_pred matrix chol(z_pred_cov)

            F77_NAME(dcopy)(&nnr, Vz, &incOne, cholVz, &incOne);

            for(k = 0; k < r; k++){
                F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
                mkLT(&cholVz[nn * k], n);
                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, &cholVz[nn * k], &n, &Cz[nn_pred * k], &n FCONE FCONE FCONE FCONE);                                  // Cz = cholinv(Vz)*Cz
                F77_NAME(dgemm)(ytran, ntran, &n_pred, &n_pred, &n, &one, &Cz[nn_pred * k], &n, &Cz[nn_pred * k], &n, &zero, &z_pred_cov[n_predn_pred * k], &n_pred FCONE FCONE);   // z_pred_cov = t(Cz)*inv(Vz)*Cz
                F77_NAME(daxpy)(&n_predn_pred, &negOne, &Vz_new[n_predn_pred * k], &incOne, &z_pred_cov[n_predn_pred * k], &incOne);
                F77_NAME(dscal)(&n_predn_pred, &negOne, &z_pred_cov[n_predn_pred * k], &incOne);
                F77_NAME(dpotrf)(lower, &n_pred, &z_pred_cov[n_predn_pred * k], &n_pred, &info FCONE); if(info != 0){perror("c++ error: z_pred_cov dpotrf failed\n");}
                mkLT(&z_pred_cov[n_predn_pred * k], n_pred);
            }

        }

        for(s = 0; s < nSamples; s++){

            F77_NAME(dcopy)(&p, &betaSamps[s * p], &incOne, beta_s, &incOne);
            F77_NAME(dcopy)(&nr, &zSamps[s * nr], &incOne, z_s, &incOne);

            if(processType == "independent.shared"){

                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &r, &one, cholVz, &n, z_s, &n FCONE FCONE FCONE FCONE);
                F77_NAME(dgemm)(ytran, ntran, &n_pred, &r, &n, &one, Cz, &n, z_s, &n, &zero, z_pred_mu, &n_pred FCONE FCONE);
                for(k = 0; k < r; k++){
                    for(i = 0; i < n_pred; i++){
                        tmp_n_predr[k * n_pred + i] = rnorm(0.0, sqrt(zScaleSamps[s]));
                    }
                }
                F77_NAME(dgemm)(ntran, ntran, &n_pred, &r, &n_pred, &one, z_pred_cov, &n_pred, tmp_n_predr, &n_pred, &zero, z_pred_s, &n_pred FCONE FCONE);
                F77_NAME(daxpy)(&n_predr, &one, z_pred_mu, &incOne, z_pred_s, &incOne);
                F77_NAME(dcopy)(&n_predr, &z_pred_s[0], &incOne, &REAL(samples_predz_r)[s*n_predr], &incOne);

            }else if(processType == "independent"){

                F77_NAME(dcopy)(&r, &zScaleSamps[s * r], &incOne, zScale_s, &incOne);
                for(k = 0; k < r; k++){
                    F77_NAME(dtrsv)(lower, ntran, nunit, &n, &cholVz[nn * k], &n, &z_s[n * k], &incOne FCONE FCONE FCONE);
                    F77_NAME(dgemv)(ytran, &n, &n_pred, &one, &Cz[nn_pred * k], &n, &z_s[n * k], &incOne, &zero, &z_pred_mu[n_pred * k], &incOne FCONE);
                    for(i = 0; i < n_pred; i++){
                        tmp_n_predr[k * n_pred + i] = rnorm(0.0, sqrt(zScale_s[k]));
                    }
                    F77_NAME(dgemv)(ntran, &n_pred, &n_pred, &one, &z_pred_cov[n_predn_pred * k], &n_pred, &tmp_n_predr[n_pred * k], &incOne, &zero, &z_pred_s[n_pred * k], &incOne FCONE);
                    F77_NAME(daxpy)(&n_pred, &one, &z_pred_mu[n_pred * k], &incOne, &z_pred_s[n_pred * k], &incOne);
                }
                F77_NAME(dcopy)(&n_predr, &z_pred_s[0], &incOne, &REAL(samples_predz_r)[s*n_predr], &incOne);

            }else if(processType == "multivariate"){

                F77_NAME(dcopy)(&rr, &zScaleSamps[s * rr], &incOne, zScale_s, &incOne);
                F77_NAME(dpotrf)(lower, &r, zScale_s, &r, &info FCONE); if(info != 0){perror("c++ error: zScale_s dpotrf failed\n");}
                mkLT(zScale_s, r);
                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &r, &one, cholVz, &n, z_s, &n FCONE FCONE FCONE FCONE);
                F77_NAME(dgemm)(ytran, ntran, &n_pred, &r, &n, &one, Cz, &n, z_s, &n, &zero, z_pred_mu, &n_pred FCONE FCONE);
                for(k = 0; k < r; k++){
                    for(i = 0; i < n_pred; i++){
                        z_pred_s[k * n_pred + i] = rnorm(0.0, 1.0);
                    }
                }
                F77_NAME(dgemm)(ntran, ntran, &n_pred, &r, &n_pred, &one, z_pred_cov, &n_pred, z_pred_s, &n_pred, &zero, tmp_n_predr, &n_pred FCONE FCONE);
                F77_NAME(dgemm)(ntran, ytran, &n_pred, &r, &r, &one, tmp_n_predr, &n_pred, zScale_s, &r, &zero, z_pred_s, &n_pred FCONE FCONE);
                F77_NAME(daxpy)(&n_predr, &one, z_pred_mu, &incOne, z_pred_s, &incOne);
                F77_NAME(dcopy)(&n_predr, &z_pred_s[0], &incOne, &REAL(samples_predz_r)[s*n_predr], &incOne);

            }

            // Find canonical parameter (X*beta + X_tilde*z_tilde)
            lmulm_XTilde_VC(ntran, n_pred, r, 1, XTilde_new, z_pred_s, tmp_n_predr);
            F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, beta_s, &incOne, &one, tmp_n_predr, &incOne FCONE);

            // Sample from predictive distribution
            for(i = 0; i < n_pred; i++){
                if(family == family_poisson){
                    dtemp1 = exp(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rpois(dtemp1);
                }else if(family == family_binary){
                    dtemp1 = inverse_logit(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rbinom(1, dtemp1);
                }else if(family == family_binomial){
                    dtemp1 = inverse_logit(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rbinom(nBinom_new[i], dtemp1);
                }
            }

        }
        // End of joint prediction set-up

    }else{

        // Implement point-wise prediction

        double *z_pred_cov = NULL;  // define NULL pointer for z_pred_cov

        // first find the Schur complement RTilde - t(C)*inv(R)*C
        // here, RTilde is always 1.0 (diagonal of correlation matrix)
        if(processType == "independent.shared" || processType == "multivariate"){

            cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);                              // nxn matrix chol(Vz)
            F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
            F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
            mkLT(cholVz, n);

            z_pred_cov = (double *) R_alloc(n_pred, sizeof(double)); zeros(z_pred_cov, n_pred);  // n_predx1 vector z_pred_cov

            F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);          // Cz = cholinv(Vz)*Cz
            for(i = 0; i < n_pred; i++){
                z_pred_cov[i] = 1.0 - F77_CALL(ddot)(&n, &Cz[i * n], &incOne, &Cz[i * n], &incOne);
            }

        }else if(processType == "independent"){

            cholVz = (double *) R_alloc(nnr, sizeof(double)); zeros(cholVz, nnr);                             // r nxn matrices chol(Vz)
            z_pred_cov = (double *) R_alloc(n_predr, sizeof(double)); zeros(z_pred_cov, n_predr);             // n_predxr vector z_pred_cov
            for(k = 0; k < r; k++){
                F77_NAME(dcopy)(&nn, &Vz[nn * k], &incOne, &cholVz[nn * k], &incOne);
                F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
                mkLT(&cholVz[nn * k], n);

                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n_pred, &one, &cholVz[nn * k], &n, &Cz[nn_pred * k], &n FCONE FCONE FCONE FCONE);          // Cz = cholinv(Vz)*Cz
                for(i = 0; i < n_pred; i++){
                    z_pred_cov[k * n_pred + i] = 1.0 - F77_CALL(ddot)(&n, &Cz[nn_pred * k + i * n], &incOne, &Cz[nn_pred * k + i * n], &incOne);
                }
            }

        }

        for(s = 0; s < nSamples; s++){

            F77_NAME(dcopy)(&p, &betaSamps[s * p], &incOne, beta_s, &incOne);
            F77_NAME(dcopy)(&nr, &zSamps[s * nr], &incOne, z_s, &incOne);

            if(processType == "independent.shared"){

                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &r, &one, cholVz, &n, z_s, &n FCONE FCONE FCONE FCONE);
                F77_NAME(dgemm)(ytran, ntran, &n_pred, &r, &n, &one, Cz, &n, z_s, &n, &zero, z_pred_mu, &n_pred FCONE FCONE);
                for(k = 0; k < r; k++){
                    for(i = 0; i < n_pred; i++){
                        z_pred_s[k * n_pred + i] = rnorm(0.0, sqrt(zScaleSamps[s] * z_pred_cov[i]));
                    }
                }
                F77_NAME(daxpy)(&n_predr, &one, z_pred_mu, &incOne, z_pred_s, &incOne);
                F77_NAME(dcopy)(&n_predr, z_pred_s, &incOne, &REAL(samples_predz_r)[s * n_predr], &incOne);

            }else if(processType == "independent"){

                F77_NAME(dcopy)(&r, &zScaleSamps[s * r], &incOne, zScale_s, &incOne);
                for(k = 0; k < r; k++){
                    F77_NAME(dtrsv)(lower, ntran, nunit, &n, &cholVz[nn * k], &n, &z_s[n * k], &incOne FCONE FCONE FCONE);
                    F77_NAME(dgemv)(ytran, &n, &n_pred, &one, &Cz[nn_pred * k], &n, &z_s[n * k], &incOne, &zero, &z_pred_mu[n_pred * k], &incOne FCONE);
                    for(i = 0; i < n_pred; i++){
                        z_pred_s[k * n_pred + i] = rnorm(0.0, sqrt(zScale_s[k] * z_pred_cov[k * n_pred + i]));
                    }
                    F77_NAME(daxpy)(&n_pred, &one, &z_pred_mu[n_pred * k], &incOne, &z_pred_s[n_pred * k], &incOne);
                }
                F77_NAME(dcopy)(&n_predr, z_pred_s, &incOne, &REAL(samples_predz_r)[s * n_predr], &incOne);

            }else if(processType == "multivariate"){

                F77_NAME(dcopy)(&rr, &zScaleSamps[s * rr], &incOne, zScale_s, &incOne);
                F77_NAME(dpotrf)(lower, &r, zScale_s, &r, &info FCONE); if(info != 0){perror("c++ error: zScale_s dpotrf failed\n");}
                mkLT(zScale_s, r);
                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &r, &one, cholVz, &n, z_s, &n FCONE FCONE FCONE FCONE);
                F77_NAME(dgemm)(ytran, ntran, &n_pred, &r, &n, &one, Cz, &n, z_s, &n, &zero, z_pred_mu, &n_pred FCONE FCONE);
                for(k = 0; k < r; k++){
                    for(i = 0; i < n_pred; i++){
                        z_pred_s[k * n_pred + i] = rnorm(0.0, sqrt(z_pred_cov[i]));
                    }
                }
                F77_NAME(dgemm)(ntran, ytran, &n_pred, &r, &r, &one, z_pred_s, &n_pred, zScale_s, &r, &zero, tmp_n_predr, &n_pred FCONE FCONE);
                F77_NAME(daxpy)(&n_predr, &one, z_pred_mu, &incOne, tmp_n_predr, &incOne);
                F77_NAME(dcopy)(&n_predr, tmp_n_predr, &incOne, &REAL(samples_predz_r)[s * n_predr], &incOne);

            }

            // Find canonical parameter (X*beta + X_tilde*z_tilde)
            lmulm_XTilde_VC(ntran, n_pred, r, 1, XTilde_new, z_pred_s, tmp_n_predr);
            F77_NAME(dgemv)(ntran, &n_pred, &p, &one, X_new, &n_pred, beta_s, &incOne, &one, tmp_n_predr, &incOne FCONE);

            // Sample from predictive distribution
            for(i = 0; i < n_pred; i++){
                if(family == family_poisson){
                    dtemp1 = exp(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rpois(dtemp1);
                }else if(family == family_binary){
                    dtemp1 = inverse_logit(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rbinom(1, dtemp1);
                }else if(family == family_binomial){
                    dtemp1 = inverse_logit(tmp_n_predr[i]);
                    REAL(samples_predmu_r)[s * n_pred + i] = dtemp1;
                    REAL(samples_predy_r)[s * n_pred + i] = rbinom(nBinom_new[i], dtemp1);
                }
            }

        }

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