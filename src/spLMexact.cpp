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

  SEXP spLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                 SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                 SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                 SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    // int i, j, k, info, nProtect = 0;
    int i, j, s, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
    // char const *rside = "R";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
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
    double dtemp = 0, dtemp2 = 0;

    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
    double *thetasp = (double *) R_alloc(2, sizeof(double));

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);
    double *tmp_np = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np, np);
    double *tmp_nn = (double *) R_alloc(nn, sizeof(double)); zeros(tmp_nn, nn);
    double *tmp_nn2 = (double *) R_alloc(nn, sizeof(double)); zeros(tmp_nn2, nn);

    double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p2, p);
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp2, pp);

    //construct covariance matrix
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull(coordsD, n, thetasp, corfn, Vz);
    F77_NAME(dcopy)(&nn, Vz, &incOne, tmp_nn, &incOne);
    for(i = 0; i < n; i++){
      tmp_nn[i*n + i] += deltasq;
    }

    // find sse
    F77_NAME(dpotrf)(lower, &n, tmp_nn, &n, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}

    F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);
    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, tmp_nn, &n, tmp_n, &incOne FCONE FCONE FCONE);

    dtemp = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne), 2);
    sse += dtemp;

    F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, tmp_nn, &n, tmp_np, &n FCONE FCONE FCONE FCONE);
    F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);

    F77_NAME(dcopy)(&pp, betaV, &incOne, tmp_pp, &incOne);
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
    F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
    F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, betaMu, &incOne, &zero, tmp_p2, &incOne FCONE);

    dtemp = F77_CALL(ddot)(&p, betaMu, &incOne, tmp_p2, &incOne);
    sse += dtemp;

    F77_NAME(daxpy)(&p, &one, tmp_p2, &incOne, tmp_p1, &incOne);
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp2, &p FCONE FCONE);
    F77_NAME(daxpy)(&pp, &one, tmp_pp2, &incOne, tmp_pp, &incOne);

    F77_NAME(dcopy)(&p, tmp_p1, &incOne, tmp_p2, &incOne);
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p2, &incOne FCONE FCONE FCONE);

    dtemp = pow(F77_NAME(dnrm2)(&p, tmp_p2, &incOne), 2);
    sse -= dtemp;

    // set-up for sampling spatial random effects
    F77_NAME(dcopy)(&nn, Vz, &incOne, tmp_nn2, &incOne);
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &n, &one, tmp_nn, &n, tmp_nn2, &n FCONE FCONE FCONE FCONE);
    F77_NAME(dgemm)(ytran, ntran, &n, &n, &n, &one, tmp_nn2, &n, tmp_nn2, &n, &zero, tmp_nn, &n FCONE FCONE);
    F77_NAME(dcopy)(&nn, Vz, &incOne, tmp_nn2, &incOne);
    F77_NAME(daxpy)(&nn, &negOne, tmp_nn, &incOne, tmp_nn2, &incOne);
    F77_NAME(dpotrf)(lower, &n, tmp_nn2, &n, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
    mkLT(tmp_nn2, n);

    // posterior parameters of sigmaSq
    sigmaSqIGaPost += sigmaSqIGa;
    sigmaSqIGaPost += 0.5 * n;

    sigmaSqIGbPost += sigmaSqIGb;
    sigmaSqIGbPost += 0.5 * sse;

    // posterior samples of sigma-sq and beta
    SEXP samples_sigmaSq_r = PROTECT(allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_beta_r = PROTECT(allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(allocMatrix(REALSXP, n, nSamples)); nProtect++;

    double sigmaSq = 0;
    double *beta = (double *) R_alloc(p, sizeof(double)); zeros(beta, p);
    double *z = (double *) R_alloc(n, sizeof(double)); zeros(z, n);

    GetRNGstate();

    for(s = 0; s < nSamples; s++){
      // sample sigmaSq from its marginal posterior
      dtemp = 1.0 / sigmaSqIGbPost;
      dtemp = rgamma(sigmaSqIGaPost, dtemp);
      sigmaSq = 1.0 / dtemp;
      REAL(samples_sigmaSq_r)[s] = sigmaSq;

      // sample fixed effects by composition sampling
      dtemp = sqrt(sigmaSq);
      for(j = 0; j < p; j++){
        beta[j] = rnorm(0.0, dtemp);
      }
      F77_NAME(daxpy)(&p, &one, tmp_p2, &incOne, beta, &incOne);
      F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, beta, &incOne FCONE FCONE FCONE);

      for(j = 0; j < p; j++){
        REAL(samples_beta_r)[s*p + j] = beta[j];
      }

      // sample spatial effects by composition sampling
      F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &one, tmp_n, &incOne FCONE);
      dtemp2 = 1.0 / deltasq;
      F77_NAME(dscal)(&n, &dtemp2, tmp_n, &incOne);

      for(i = 0; i < n; i++){
        z[i] = rnorm(0.0, dtemp);
      }

      F77_NAME(dgemv)(ytran, &n, &n, &one, tmp_nn2, &n, tmp_n, &incOne, &one, z, &incOne FCONE);
      F77_NAME(dgemv)(ntran, &n, &n, &one, tmp_nn2, &n, z, &incOne, &zero, tmp_n, &incOne FCONE);

      for(i = 0; i < n; i++){
        REAL(samples_z_r)[s*n + i] = tmp_n[i];
      }

    }

    PutRNGstate();

    // make return object for posterior samples of sigma-sq and beta
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    result_r = PROTECT(allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta"));

    // samples of sigma-sq
    SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("sigmaSq"));

    // samples of z
    SET_VECTOR_ELT(result_r, 2, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z"));

    namesgets(result_r, resultName_r);

    // SEXP samples_beta_r = PROTECT(allocMatrix(REALSXP, nSamples, p)); nProtect++;
    // SEXP tmp_n_r = PROTECT(allocVector(REALSXP, n)); nProtect++;
    // SEXP tmp_nn_r = PROTECT(allocMatrix(REALSXP, n, n)); nProtect++;

    // for (i = 0; i < n; i++) {
    //   for (j = 0; j < n; j++) {
    //     REAL(tmp_nn_r)[i*n + j] = tmp_nn2[i*n + j];
    //   }
    // }

    // for(i = 0; i < n; i++){
    //   REAL(tmp_n_r)[i] = tmp_n[i];
    // }

    // REAL(result_r)[0] = sigmaSqIGaPost;
    // REAL(result_r)[1] = sigmaSqIGbPost;

    // REAL(result_r)[0] = sse;

    UNPROTECT(nProtect);

    return result_r;

  }

}
