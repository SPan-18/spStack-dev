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

  SEXP recoverScale_stvcGLM(SEXP n_r, SEXP p_r, SEXP r_r, SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                            SEXP betaMu_r, SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP iwScale_r, SEXP processType_r,
                            SEXP phi_s_r, SEXP phi_t_r, SEXP nSamples_r, SEXP betaSamps_r, SEXP zSamps_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, k, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
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
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int r = INTEGER(r_r)[0];
    int rr = r * r;
    int nr = n * r;
    int nnr = nn * r;
    double *zSamps = REAL(zSamps_r);
    double *betaSamps = REAL(betaSamps_r);

    double *coords_sp = REAL(sp_coords_r);
    double *coords_tm = REAL(time_coords_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    double *betaMu = (double *) R_alloc(p, sizeof(double)); zeros(betaMu, p);
    F77_NAME(dcopy)(&p, REAL(betaMu_r), &incOne, betaMu, &incOne);
    double *betaV = (double *) R_alloc(pp, sizeof(double)); zeros(betaV, pp);
    F77_NAME(dcopy)(&pp, REAL(betaV_r), &incOne, betaV, &incOne);

    double nu_beta = REAL(nu_beta_r)[0];
    double nu_z = REAL(nu_z_r)[0];

    double *iwScale = (double *) R_alloc(rr, sizeof(double)); zeros(iwScale, rr);
    F77_NAME(dcopy)(&rr, REAL(iwScale_r), &incOne, iwScale, &incOne);

    // spatial-temporal process parameters: create spatial-temporal covariance matrices
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    double *Vz = NULL;

    if(corfn == "gneiting-decay"){

        if(processType == "independent.shared" || processType == "multivariate"){

            phi_s_vec[0] = REAL(phi_s_r)[0];
            phi_t_vec[0] = REAL(phi_t_r)[0];
            thetaspt[0] = phi_s_vec[0];
            thetaspt[1] = phi_t_vec[0];

            Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
            sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, Vz);

        }else if(processType == "independent"){

            F77_NAME(dcopy)(&r, REAL(phi_s_r), &incOne, phi_s_vec, &incOne);
            F77_NAME(dcopy)(&r, REAL(phi_t_r), &incOne, phi_t_vec, &incOne);

            Vz = (double *) R_alloc(nnr, sizeof(double)); zeros(Vz, nnr);

            // find r-many correlation matrices, stacked into a rn^2-dim vector
            for(k = 0; k < r; k++){
                thetaspt[0] = phi_s_vec[k];
                thetaspt[1] = phi_t_vec[k];
                sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, &Vz[nn * k]);
            }

        }
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/

    double *cholVz = NULL;               // define NULL pointer for chol(Vz)

    // Find Cholesky of Vz
    if(processType == "independent.shared"){

        cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);            // nxn matrix chol(Vz)
        F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
        mkLT(cholVz, n);

    }else if(processType == "independent"){

        cholVz = (double *) R_alloc(nnr, sizeof(double)); zeros(cholVz, nnr);          // r nxn matrices chol(Vz)
        F77_NAME(dcopy)(&nnr, Vz, &incOne, cholVz, &incOne);
        for(k = 0; k < r; k++){
            F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
            mkLT(&cholVz[nn * k], n);
        }

    }else if(processType == "multivariate"){

        cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);            // nxn matrix chol(Vz)
        F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
        mkLT(cholVz, n);

    }

    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                                    // Cholesky of Vbeta
    F77_NAME(dcopy)(&pp, betaV, &incOne, Lbeta, &incOne);                                                        // Lbeta = Vbeta
    F77_NAME(dpotrf)(lower, &p, Lbeta, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}    // Lbeta = chol(Vbeta)

    /*****************************************
     Set-up posterior sampling
     *****************************************/
    // posterior samples of sigma-sq and beta
    SEXP samples_betaScale_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_zScale_r = R_NilValue;
    if(processType == "independent.shared"){
        samples_zScale_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    }else if(processType == "independent"){
        samples_zScale_r = PROTECT(Rf_allocMatrix(REALSXP, r, nSamples)); nProtect++;
    }else if(processType == "multivariate"){
        samples_zScale_r = PROTECT(Rf_allocMatrix(REALSXP, rr, nSamples)); nProtect++;
    }

    // temmporary variables for posterior recovery of scale parameters
    double QBeta = 0.0, Qz = 0.0, IGa = 0.0, IGb = 0.0;
    double *beta_s = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta_s, p);
    double *z_s = (double *) R_chk_calloc(nr, sizeof(double)); zeros(z_s, nr);
    double *Qz_rr = (double *) R_chk_calloc(rr, sizeof(double)); zeros(Qz_rr, rr);
    double *samp_Sigma = (double *) R_chk_calloc(rr, sizeof(double)); zeros(samp_Sigma, rr);
    double *tmp_rr = (double *) R_chk_calloc(rr, sizeof(double)); zeros(tmp_rr, rr);

    // recover posterior samples of scale parameter of beta
    for(i = 0; i < nSamples; i++){
        F77_NAME(dcopy)(&p, &betaSamps[i * p], &incOne, beta_s, &incOne);
        F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, beta_s, &incOne);
        F77_NAME(dtrsv)(lower, ntran, nUnit, &p, Lbeta, &p, beta_s, &incOne FCONE FCONE FCONE);
        QBeta = F77_NAME(ddot)(&p, beta_s, &incOne, beta_s, &incOne);
        IGa = 0.5 * (nu_beta + p);
        IGb = 0.5 * (nu_beta + QBeta);
        REAL(samples_betaScale_r)[i] = 1.0 / rgamma(IGa, 1.0 / IGb);
    }

    // recover posterior samples of scale parameter of z
    for(i = 0; i < nSamples; i++){

        if(processType == "independent.shared"){

            F77_NAME(dcopy)(&nr, &zSamps[i * nr], &incOne, z_s, &incOne);
            Qz = 0.0;
            for(j = 0; j < r; j++){
                F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVz, &n, &z_s[n * j], &incOne FCONE FCONE FCONE);
                Qz += F77_NAME(ddot)(&n, &z_s[n * j], &incOne, &z_s[n * j], &incOne);
            }
            IGa = 0.5 * (nu_z + nr);
            IGb = 0.5 * (nu_z + Qz);
            REAL(samples_zScale_r)[i] = 1.0 / rgamma(IGa, 1.0 / IGb);

        }else if(processType == "independent"){

            F77_NAME(dcopy)(&nr, &zSamps[i * nr], &incOne, z_s, &incOne);
            for(j = 0; j < r; j++){
                F77_NAME(dtrsv)(lower, ntran, nUnit, &n, &cholVz[nn * j], &n, &z_s[n * j], &incOne FCONE FCONE FCONE);
                Qz = F77_NAME(ddot)(&n, &z_s[n * j], &incOne, &z_s[n * j], &incOne);
                IGa = 0.5 * (nu_z + nr);
                IGb = 0.5 * (nu_z + Qz);
                REAL(samples_zScale_r)[i * r + j] = 1.0 / rgamma(IGa, 1.0 / IGb);
            }

        }else if(processType == "multivariate"){

            F77_NAME(dcopy)(&nr, &zSamps[i * nr], &incOne, z_s, &incOne);
            F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &r, &one, cholVz, &n, z_s, &n FCONE FCONE FCONE FCONE);                // z_s = cholinv(R)*Z
            F77_NAME(dgemm)(ytran, ntran, &r, &r, &n, &one, z_s, &n, z_s, &n, &zero, Qz_rr, &r FCONE FCONE);                       // Qz = Z'inv(R)Z
            F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, Qz_rr, &incOne);                                                          // Qz = Z'inv(R)Z + iwScale
            F77_NAME(dpotrf)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");} // chol(Qz)
            F77_NAME(dpotri)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotri failed\n");} // inv(Qz)
            F77_NAME(dpotrf)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");} // chol(inv(Qz))
            mkLT(Qz_rr, r);
            rInvWishart(r, nu_z + n + 2*r, Qz_rr, samp_Sigma, tmp_rr);
            F77_NAME(dcopy)(&rr, samp_Sigma, &incOne, &REAL(samples_zScale_r)[i * rr], &incOne);

        }
    }

    R_chk_free(beta_s);
    R_chk_free(z_s);
    R_chk_free(Qz_rr);
    R_chk_free(samp_Sigma);
    R_chk_free(tmp_rr);

    // make return object
    SEXP result_r, resultName_r;

    // make return object for posterior samples
    int nResultListObjs = 2;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of scale parameter of beta
    SET_VECTOR_ELT(result_r, 0, samples_betaScale_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("sigmasq.beta"));

    // samples of scale parameter of z
    SET_VECTOR_ELT(result_r, 1, samples_zScale_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("z.scale"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

  }  // End of recoverScale_stvcGLM

}