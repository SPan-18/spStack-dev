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

    SEXP recoverScale_spGLM(SEXP n_r, SEXP p_r, SEXP coordsD_r, SEXP corfn_r,
                            SEXP betaMu_r, SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r,
                            SEXP phi_r, SEXP nu_r, SEXP nSamples_r, SEXP betaSamps_r, SEXP zSamps_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    const double negOne = -1.0;
    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    double *zSamps = REAL(zSamps_r);
    double *betaSamps = REAL(betaSamps_r);

    double *coordsD = REAL(coordsD_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    double *betaMu = (double *) R_alloc(p, sizeof(double)); zeros(betaMu, p);
    F77_NAME(dcopy)(&p, REAL(betaMu_r), &incOne, betaMu, &incOne);
    double *betaV = (double *) R_alloc(pp, sizeof(double)); zeros(betaV, pp);
    F77_NAME(dcopy)(&pp, REAL(betaV_r), &incOne, betaV, &incOne);

    double nu_beta = REAL(nu_beta_r)[0];
    double nu_z = REAL(nu_z_r)[0];

    // spatial process parameters
    double phi = REAL(phi_r)[0];

    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];

    // memory allocations
    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);                       // correlation matrix
    double *thetasp = (double *) R_alloc(2, sizeof(double));                                  // spatial process parameters
    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);               // Cholesky of Vz
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                 // Cholesky of Vbeta

    //construct covariance matrix (full)
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull(coordsD, n, thetasp, corfn, Vz);

    // Find Cholesky of Vz
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}

    F77_NAME(dcopy)(&pp, betaV, &incOne, Lbeta, &incOne);                                                            // Lbeta = Vbeta
    F77_NAME(dpotrf)(lower, &p, Lbeta, &p, &info FCONE); if(info != 0){perror("c++ error: VBeta dpotrf failed\n");}  // Lbeta = chol(Vbeta)

    /*****************************************
     Set-up posterior sampling
     *****************************************/
    // posterior samples of sigma-sq and beta
    SEXP samples_betaScale_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_zScale_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;

    // temmporary variables for posterior recovery of scale parameters
    double QBeta = 0.0, Qz = 0.0, IGa = 0.0, IGb = 0.0;
    double *beta_s = (double *) R_alloc(p, sizeof(double)); zeros(beta_s, p);
    double *z_s = (double *) R_alloc(n, sizeof(double)); zeros(z_s, n);

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
        F77_NAME(dcopy)(&n, &zSamps[i * n], &incOne, z_s, &incOne);
        F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVz, &n, z_s, &incOne FCONE FCONE FCONE);
        Qz = F77_NAME(ddot)(&n, z_s, &incOne, z_s, &incOne);
        IGa = 0.5 * (nu_z + n);
        IGb = 0.5 * (nu_z + Qz);
        REAL(samples_zScale_r)[i] = 1.0 / rgamma(IGa, 1.0 / IGb);
    }

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
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmasq.z"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

    }

}