#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP spLM_fixed(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                  SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                  SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                  SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    // int i, j, k, info, nProtect = 0;
    int i, j, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
    // char const *rside = "R";
    char const *lside = "L";
    const double one = 1.0;
    // const double negOne = -1.0;
    const double zero = 0.0;
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

    double *coordsD = REAL(coordsD_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r, 0));
    double *betaMu = NULL;
    double *betaV = NULL;

    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);

      betaV = (double *) R_alloc(pp, sizeof(double));
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaV, &incOne);
    }

    double sigmaSqIGa = REAL(sigmaSqIG_r)[0];
    double sigmaSqIGb = REAL(sigmaSqIG_r)[1];

    double deltasq = REAL(deltasq_r)[0];
    double phi = REAL(phi_r)[0];

    double nu = 0;
    if(corfn == "matern"){
      nu = REAL(nu_r)[0];
    }

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];

    // print set-up if verbose TRUE
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept).\n\n", p);
      Rprintf("Using the %s spatial correlation function.\n\n", corfn.c_str());

      Rprintf("Priors:\n");

      if(betaPrior == "flat"){
        Rprintf("\tbeta flat.\n");
      }else{
        Rprintf("\tbeta normal:\n");
        Rprintf("\tmu:"); printVec(betaMu, p);
        Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
        Rprintf("\n");
      }

      Rprintf("\tsigma.sq IG hyperpriors shape = %.5f and scale = %.5f\n\n", sigmaSqIGa, sigmaSqIGb);

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.5f, and, nu = %.5f\n", phi, nu);
      }else{
        Rprintf("\tphi = %.5f\n", phi);
      }
      Rprintf("\tNoise-to-spatial variance ratio = %.5f\n\n", deltasq);

      Rprintf("Number of posterior samples = %i.\n", nSamples);

    }

    /*****************************************
     Set-up posterior sample vector/matrices etc.
     *****************************************/
    double sigmaSqIGaPost = 0, sigmaSqIGbPost = 0;
    double sse = 0;

    double *V = (double *) R_alloc(nn, sizeof(double)); zeros(V, nn);
    double *Vy = (double *) R_alloc(nn, sizeof(double)); zeros(Vy, nn);
    double *thetasp = (double *) R_alloc(2, sizeof(double));

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);
    double *tmp_np = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np, np);
    double *tmp_nn = (double *) R_alloc(nn, sizeof(double)); zeros(tmp_nn, nn);

    double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p2, p);
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp2, pp);

    //construct covariance matrix
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorLT(coordsD, n, thetasp, corfn, V);
    F77_NAME(dcopy)(&nn, V, &incOne, Vy, &incOne);
    for(i = 0; i < n; i++){
      Vy[i*n + i] = V[i*n + i] + deltasq;
    }

    F77_NAME(dcopy)(&nn, Vy, &incOne, tmp_nn, &incOne);
    F77_NAME(dpotrf)(lower, &n, tmp_nn, &n, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}

    F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);
    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, tmp_nn, &n, tmp_n, &incOne FCONE FCONE FCONE);

    sse += pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne), 2);

    F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, tmp_nn, &n, tmp_np, &n FCONE FCONE FCONE FCONE);
    F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);

    F77_NAME(dcopy)(&pp, betaV, &incOne, tmp_pp, &incOne);
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
    F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
    F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, betaMu, &incOne, &zero, tmp_p2, &incOne FCONE);

    double ssetemp = F77_CALL(ddot)(&p, betaMu, &incOne, tmp_p2, &incOne);
    sse += ssetemp;

    F77_NAME(daxpy)(&p, &one, tmp_p2, &incOne, tmp_p1, &incOne);
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp2, &p FCONE FCONE);

    for(i = 0; i < p; i++){
      for(j = 0; j < p; j++){
        tmp_pp[i*p + j] += tmp_pp2[i*p + j];
      }
    }

    F77_NAME(dcopy)(&p, tmp_p1, &incOne, tmp_p2, &incOne);
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p2, &incOne FCONE FCONE FCONE);
    // mysolveLT(tmp_pp, tmp_p2, p);

    sse -= pow(F77_NAME(dnrm2)(&p, tmp_p2, &incOne), 2);

    // posterior parameters of sigmaSq
    sigmaSqIGaPost += sigmaSqIGa;
    sigmaSqIGaPost += 0.5 * n;

    sigmaSqIGbPost += sigmaSqIGb;
    sigmaSqIGbPost += 0.5 * sse;

    // posterior samples of sigma-sq
    double *sigmaSq = (double *) R_alloc(nSamples, sizeof(double));
    for(i = 0; i < nSamples; i++){
      sigmaSq[i] = rgamma(sigmaSqIGaPost, 1.0 / sigmaSqIGbPost);
      sigmaSq[i] = 1.0 / sigmaSq[i];
    }

    //return posterior samples of sigma-sq
    SEXP samples_sigmaSq_r = PROTECT(allocVector(REALSXP, nSamples)); nProtect++;
    // SEXP result_r = PROTECT(allocVector(REALSXP, 2)); nProtect++;
    // SEXP result_r = PROTECT(allocMatrix(REALSXP, p, p));

    // for (i = 0; i < p; i++) {
    //   for (j = 0; j < p; j++) {
    //     REAL(result_r)[i*p + j] = tmp_pp[i*p + j];
    //   }
    // }

    for(i = 0; i < nSamples; i++){
      REAL(samples_sigmaSq_r)[i] = sigmaSq[i];
    }

    // REAL(result_r)[0] = sigmaSqIGaPost;
    // REAL(result_r)[1] = sigmaSqIGbPost;

    // REAL(result_r)[0] = sse;

    UNPROTECT(nProtect);

    return samples_sigmaSq_r;

  }

}
