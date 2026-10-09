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

    // supported spatial-temporal process models
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    if(processType != "independent.shared" && processType != "independent" && processType != "multivariate"){
      Rf_error("c++ error: process.type must be one of 'independent', 'independent.shared' or 'multivariate'.");
    }
    if(corfn != "gneiting-decay"){
      Rf_error("c++ error: cor.fn must be 'gneiting-decay'.");
    }
    int nCov = (processType == "independent") ? r : 1;                   // number of distinct correlation functions

    // spatial-temporal process parameters and chol(R_k), built in place
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double thetaspt[2] = {0.0, 0.0};
    F77_NAME(dcopy)(&nCov, REAL(phi_s_r), &incOne, phi_s_vec, &incOne);
    F77_NAME(dcopy)(&nCov, REAL(phi_t_r), &incOne, phi_t_vec, &incOne);
    double *cholVz = (double *) R_alloc((size_t) nn * nCov, sizeof(double)); zeros(cholVz, nn * nCov);
    for(k = 0; k < nCov; k++){
        thetaspt[0] = phi_s_vec[k];
        thetaspt[1] = phi_t_vec[k];
        sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, &cholVz[nn * k]);
        F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE);
        if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i).\n", info);}
        mkLT(&cholVz[nn * k], n);
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];

    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                                    // Cholesky of Vbeta
    F77_NAME(dcopy)(&pp, betaV, &incOne, Lbeta, &incOne);                                                        // Lbeta = Vbeta
    F77_NAME(dpotrf)(lower, &p, Lbeta, &p, &info FCONE);                                                         // Lbeta = chol(Vbeta)
    if(info != 0){Rf_error("c++ error: Cholesky factorization of the prior covariance of beta failed (info = %i).\n", info);}

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
    }else{
        samples_zScale_r = PROTECT(Rf_allocMatrix(REALSXP, rr, nSamples)); nProtect++;
    }

    // temmporary variables for posterior recovery of scale parameters
    double QBeta = 0.0, Qz = 0.0, IGa = 0.0, IGb = 0.0;
    double *beta_s = (double *) R_alloc(p, sizeof(double)); zeros(beta_s, p);
    double *Qz_rr = (double *) R_alloc(rr, sizeof(double)); zeros(Qz_rr, rr);
    double *samp_Sigma = (double *) R_alloc(rr, sizeof(double)); zeros(samp_Sigma, rr);
    double *tmp_rr = (double *) R_alloc(rr, sizeof(double)); zeros(tmp_rr, rr);

    // the solves cholinv(R_k)*z_k are done for blocks of nBlockMax samples with level-3 BLAS; the inverse-gamma
    // and inverse-Wishart draws are made in the original order
    const int nBlockMax = 64;
    int nBlock = 0, nrBlock = 0, rBlock = 0, b = 0, s = 0;
    double *zBlock = (double *) R_alloc((size_t) nr * nBlockMax, sizeof(double));            // nr x nBlock work matrix
    double *Zb = NULL;

    GetRNGstate();

    // recover posterior samples of scale parameter of beta
    for(i = 0; i < nSamples; i++){
        F77_NAME(dcopy)(&p, &betaSamps[(R_xlen_t) i * p], &incOne, beta_s, &incOne);
        F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, beta_s, &incOne);
        F77_NAME(dtrsv)(lower, ntran, nUnit, &p, Lbeta, &p, beta_s, &incOne FCONE FCONE FCONE);
        QBeta = F77_NAME(ddot)(&p, beta_s, &incOne, beta_s, &incOne);
        IGa = 0.5 * (nu_beta + p);
        IGb = 0.5 * (nu_beta + QBeta);
        REAL(samples_betaScale_r)[i] = 1.0 / rgamma(IGa, 1.0 / IGb);
    }

    // recover posterior samples of scale parameter of z
    for(s = 0; s < nSamples; s += nBlockMax){

        nBlock = std::min(nBlockMax, nSamples - s);
        nrBlock = nr * nBlock;
        rBlock = r * nBlock;
        F77_NAME(dcopy)(&nrBlock, &zSamps[(R_xlen_t) s * nr], &incOne, zBlock, &incOne);
        if(nCov == 1){
            F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &rBlock, &one, cholVz, &n, zBlock, &n FCONE FCONE FCONE FCONE);   // cholinv(R)*z_j, all j
        }else{
            for(j = 0; j < r; j++){
                F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &nBlock, &one, &cholVz[nn * j], &n, &zBlock[n * j], &nr FCONE FCONE FCONE FCONE);
            }
        }

        for(b = 0; b < nBlock; b++){

            i = s + b;
            Zb = &zBlock[nr * b];

            if(processType == "independent.shared"){

                Qz = 0.0;
                for(j = 0; j < r; j++){
                    Qz += F77_NAME(ddot)(&n, &Zb[n * j], &incOne, &Zb[n * j], &incOne);
                }
                IGa = 0.5 * (nu_z + nr);
                IGb = 0.5 * (nu_z + Qz);
                REAL(samples_zScale_r)[i] = 1.0 / rgamma(IGa, 1.0 / IGb);

            }else if(processType == "independent"){

                for(j = 0; j < r; j++){
                    Qz = F77_NAME(ddot)(&n, &Zb[n * j], &incOne, &Zb[n * j], &incOne);
                    IGa = 0.5 * (nu_z + n);                                    // process j has its own scale and n values
                    IGb = 0.5 * (nu_z + Qz);
                    REAL(samples_zScale_r)[(R_xlen_t) i * r + j] = 1.0 / rgamma(IGa, 1.0 / IGb);
                }

            }else{

                F77_NAME(dgemm)(ytran, ntran, &r, &r, &n, &one, Zb, &n, Zb, &n, &zero, Qz_rr, &r FCONE FCONE);                        // Qz = Z'inv(R)Z
                F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, Qz_rr, &incOne);                                                          // Qz = Z'inv(R)Z + iwScale
                F77_NAME(dpotrf)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");} // chol(Qz)
                F77_NAME(dpotri)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotri failed\n");} // inv(Qz)
                F77_NAME(dpotrf)(lower, &r, Qz_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");} // chol(inv(Qz))
                mkLT(Qz_rr, r);
                rInvWishart(r, nu_z + n + 2*r, Qz_rr, samp_Sigma, tmp_rr);
                F77_NAME(dcopy)(&rr, samp_Sigma, &incOne, &REAL(samples_zScale_r)[(R_xlen_t) i * rr], &incOne);

            }
        }
    }

    PutRNGstate();

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