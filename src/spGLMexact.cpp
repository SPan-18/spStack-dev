#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
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

  SEXP spGLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r,
                  SEXP coordsD_r, SEXP corfn_r, SEXP betaV_r, SEXP nu_beta_r,
                  SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP phi_r, SEXP nu_r,
                  SEXP epsilon_r, SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    // int i, j, s, info, nProtect = 0;
    int nProtect = 0;
    // char const *lower = "L";
    // char const *nUnit = "N";
    // char const *ntran = "N";
    // char const *ytran = "T";
    // char const *lside = "L";
    // const double one = 1.0;
    // const double negOne = -1.0;
    // const double zero = 0.0;
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
    int np = n * p;

    std::string family = CHAR(STRING_ELT(family_r, 0));

    double *coordsD = REAL(coordsD_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    double *betaMu = (double *) R_alloc(p, sizeof(double)); zeros(betaMu, p);
    double *betaV = (double *) R_alloc(pp, sizeof(double)); zeros(betaV, pp);
    F77_NAME(dcopy)(&pp, REAL(betaV_r), &incOne, betaV, &incOne);

    double nu_beta = REAL(nu_beta_r)[0];
    double nu_z = REAL(nu_z_r)[0];
    double sigmaSq_xi = REAL(sigmaSq_xi_r)[0];

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

      Rprintf("\tbeta normal:\n");
      Rprintf("\tmu:"); printVec(betaMu, p);
      Rprintf("\tcov:"); printMtrx(betaV, p, p);
      Rprintf("\n");

      Rprintf("\tsigmaSq.beta ~ IG(nu.beta/2, nu.beta/2),\tnu.beta = %.2f\n", nu_beta);
      Rprintf("\tsigmaSq.z ~ IG(nu.z/2, nu.z/2),\tnu.z = %.2f\n", nu_z);
      Rprintf("\tsigmaSq.xi = %.2f\n\n", sigmaSq_xi);

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.3f, and, nu = %.3f\n\n", phi, nu);
      }else{
        Rprintf("\tphi = %.3f\n\n", phi);
      }

      Rprintf("Number of posterior samples = %i.\n", nSamples);

    }

    SEXP result_r;
    result_r = PROTECT(Rf_allocVector(REALSXP, 1)); nProtect++;

    REAL(result_r)[0] = 0.1;

    UNPROTECT(nProtect);

    return result_r;

  } // end spGLMexact
}