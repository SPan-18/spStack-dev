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

  SEXP spGLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP nBinom_r,
                  SEXP coords_r, SEXP corfn_r, SEXP betaV_r, SEXP nu_beta_r,
                  SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP phi_r, SEXP nu_r,
                  SEXP epsilon_r, SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
    char const *lower = "L";
    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    double *Y = REAL(Y_r);
    double *nBinom = REAL(nBinom_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int np = n * p;

    std::string family = CHAR(STRING_ELT(family_r, 0));

    double *coords = REAL(coords_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    double *betaMu = (double *) R_alloc(p, sizeof(double)); zeros(betaMu, p);
    double *betaV = (double *) R_alloc(pp, sizeof(double)); zeros(betaV, pp);
    F77_NAME(dcopy)(&pp, REAL(betaV_r), &incOne, betaV, &incOne);

    double nu_beta = REAL(nu_beta_r)[0];
    double nu_z = REAL(nu_z_r)[0];
    double sigmaSq_xi = REAL(sigmaSq_xi_r)[0];
    double sigma_xi = sqrt(sigmaSq_xi);

    // spatial process parameters
    double phi = REAL(phi_r)[0];

    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    // boundary adjustment parameter
    double epsilon = REAL(epsilon_r)[0];

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];

    // print set-up if verbose TRUE
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Family = %s.\n\n", family.c_str());
      Rprintf("Number of covariates %i (including intercept).\n\n", p);
      Rprintf("Using the %s spatial correlation function.\n\n", corfn.c_str());

      Rprintf("Priors:\n");

      Rprintf("\tbeta: Gaussian\n");
      Rprintf("\tmu:"); printVec(betaMu, p);
      Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
      Rprintf("\n");

      Rprintf("\tsigmaSq.beta ~ IG(nu.beta/2, nu.beta/2)\n");
      Rprintf("\tsigmaSq.z ~ IG(nu.z/2, nu.z/2)\n");
      Rprintf("\tnu.beta = %.2f, nu.z = %.2f.\n", nu_beta, nu_z);
      Rprintf("\tsigmaSq.xi = %.2f.\n", sigmaSq_xi);
      Rprintf("\tBoundary adjustment parameter = %.2f.\n\n", epsilon);

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.2f, and, nu = %.2f.\n\n", phi, nu);
      }else{
        Rprintf("\tphi = %.2f.\n\n", phi);
      }

      Rprintf("Number of posterior samples = %i.\n", nSamples);
      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/
    double dtemp1, dtemp2, dtemp3;

    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);               // correlation matrix Vz, then its Cholesky
    double *cholVzPlusI = (double *) R_alloc(nn, sizeof(double)); zeros(cholVzPlusI, nn);     // allocate memory for n x n matrix
    double *cholSchur_n = (double *) R_alloc(nn, sizeof(double)); zeros(cholSchur_n, nn);     // allocate memory for Schur complement
    double *cholSchur_p = (double *) R_alloc(pp, sizeof(double)); zeros(cholSchur_p, pp);     // allocate memory for Schur complement
    double *D1invX = (double *) R_alloc(np, sizeof(double)); zeros(D1invX, np);               // allocate for preprocessing
    double *DinvB_pn = (double *) R_alloc(np, sizeof(double)); zeros(DinvB_pn, np);           // allocate memory for n x p matrix DinvB_np
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);           // allocate VbetaInv
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                 // Cholesky of Vbeta
    double *thetasp = (double *) R_alloc(2, sizeof(double));                                  // spatial process parameters

    //construct covariance matrix (full)
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull2(n, 2, coords, thetasp, corfn, cholVz);                                       // cholVz = Vz

    // construct unit spherical perturbation of Vz; (Vz+I)
    F77_NAME(dcopy)(&nn, cholVz, &incOne, cholVzPlusI, &incOne);
    for(i = 0; i < n; i++){
      cholVzPlusI[i*n + i] += 1.0;
    }

    // find Cholesky factor of unit spherical perturbation of Vz
    F77_NAME(dpotrf)(lower, &n, cholVzPlusI, &n, &info FCONE);
    if(info != 0){Rf_error("c++ error: Cholesky factorization of Vz + I failed (info = %i).\n", info);}

    // Find Cholesky of Vz
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE);
    if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial correlation matrix failed (info = %i); it is numerically singular, check for nearly coincident locations or a very small phi.\n", info);}

    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                           // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){Rf_error("c++ error: prior covariance of beta is not positive definite.\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);                                                           // Lbeta = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){Rf_error("c++ error: inversion of the prior covariance of beta failed.\n");} // VbetaInv = chol2inv(Vbeta)

    // Get the Schur complement of top left nxn submatrix of (HtH)
    double *tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np, np);       // temporary allocate memory for n x p matrix

    info = cholSchurGLM(X, n, p, sigmaSq_xi, VbetaInv, cholVzPlusI, tmp_np,
                        DinvB_pn, cholSchur_p, cholSchur_n, D1invX);

    R_chk_free(tmp_np);
    if(info != 0){glmPrimingError(info);}

    /*****************************************
     Set-up posterior sampling
     *****************************************/
    // posterior samples of sigma-sq and beta
    SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;
    SEXP samples_xi_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

    const char *family_poisson = "poisson";
    const char *family_binary = "binary";
    const char *family_binomial = "binomial";

    // Posterior samples are drawn in blocks of nBlockMax: within a block the random variates are drawn in
    // exactly the order of a draw-by-draw loop, directly into the output matrices, and the block is then
    // projected at once with level-3 BLAS (projGLMbatch).
    const int nBlockMax = 64;
    int nBlock = 0, bb = 0;
    int nnBlockMax = n * nBlockMax;
    int pnBlockMax = p * nBlockMax;
    double *V_eta = (double *) R_chk_calloc(nnBlockMax, sizeof(double)); zeros(V_eta, nnBlockMax);
    double *tmp_nb = (double *) R_chk_calloc(nnBlockMax, sizeof(double)); zeros(tmp_nb, nnBlockMax);
    double *tmp_pb = (double *) R_chk_calloc(pnBlockMax, sizeof(double)); zeros(tmp_pb, pnBlockMax);
    double *V_beta = NULL, *V_z = NULL, *V_xi = NULL;

    GetRNGstate();

    for(s = 0; s < nSamples; s += nBlockMax){

      nBlock = std::min(nBlockMax, nSamples - s);
      V_beta = &REAL(samples_beta_r)[(R_xlen_t) s * p];
      V_xi = &REAL(samples_xi_r)[(R_xlen_t) s * n];
      V_z = &REAL(samples_z_r)[(R_xlen_t) s * n];

      for(bb = 0; bb < nBlock; bb++){


        if(family == family_poisson){
          for(i = 0; i < n; i++){
            dtemp1 = Y[i] + epsilon;
            dtemp2 = 1.0;
            V_eta[bb*n + i] = rlogGamma(dtemp1);                              // log(Gamma(y + epsilon, 1)), underflow-safe
          }
        }

        if(family == family_binomial){
          for(i = 0; i < n; i++){
            dtemp1 = Y[i] + epsilon;
            dtemp2 = nBinom[i];
            dtemp2 += 2.0 * epsilon;
            dtemp2 -= dtemp1;
            V_eta[bb*n + i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
          }
        }

        if(family == family_binary){
          for(i = 0; i < n; i++){
            dtemp1 = Y[i] + epsilon;
            dtemp2 = nBinom[i];
            dtemp2 += 2.0 * epsilon;
            dtemp2 -= dtemp1;
            V_eta[bb*n + i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
          }
        }

        dtemp1 = 0.5 * nu_beta;
        dtemp2 = 1.0 / dtemp1;
        dtemp3 = rgamma(dtemp1, dtemp2);
        dtemp3 = 1.0 / dtemp3;
        dtemp3 = sqrt(dtemp3);
        for(j = 0; j < p; j++){
          V_beta[bb*p + j] = rnorm(0.0, dtemp3);                                                  // v_beta ~ N(0, 1)
        }

        dtemp1 = 0.5 * nu_z;
        dtemp2 = 1.0 / dtemp1;
        dtemp3 = rgamma(dtemp1, dtemp2);
        dtemp3 = 1.0 / dtemp3;
        dtemp3 = sqrt(dtemp3);
        for(i = 0; i < n; i++){
          V_xi[bb*n + i] = rnorm(0.0, sigma_xi);                                                  // v_xi ~ N(0, 1)
          V_z[bb*n + i] = rnorm(0.0, dtemp3);                                                     // v_z ~ N(0, 1)
        }

      }

      // projection step for the block
      projGLMbatch(X, n, p, nBlock, V_eta, V_xi, V_beta, V_z, cholSchur_p, cholSchur_n, sigmaSq_xi, Lbeta, cholVz, cholVzPlusI, D1invX, DinvB_pn,
                   tmp_nb, tmp_pb);

    }

    PutRNGstate();

    R_chk_free(V_eta);
    R_chk_free(tmp_nb);
    R_chk_free(tmp_pb);

    // make return object
    SEXP result_r, resultName_r;

    // make return object for posterior samples and leave-one-out predictive densities
    int nResultListObjs = 3;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

    // samples of z
    SET_VECTOR_ELT(result_r, 1, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("z"));

    // samples of z
    SET_VECTOR_ELT(result_r, 2, samples_xi_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("xi"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

  } // end spGLMexact
}