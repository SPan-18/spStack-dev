#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
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

  SEXP spLMexact2(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                  SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                  SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                  SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
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
    // int nPp = n + p;

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
    const double delta = sqrt(deltasq);
    const double deltaInv = 1.0 / delta;
    // const double deltasqInv = 1.0 / deltasq;

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
    double *Lz = (double *) R_alloc(nn, sizeof(double)); zeros(Lz, nn);
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);
    double *cholVy = (double *) R_alloc(nn, sizeof(double)); zeros(cholVy, nn);
    double *thetasp = (double *) R_alloc(2, sizeof(double));

    double *tmp_n1 = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n1, n);
    double *tmp_n2 = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n2, n);

    double *tmp_np1 = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np1, np);
    // double *tmp_np2 = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np2, np);

    double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p2, p);
    double *mu_vbeta = (double *) R_alloc(p, sizeof(double)); zeros(mu_vbeta, p);

    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);

    //construct covariance matrix
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull(coordsD, n, thetasp, corfn, Vz);
    F77_NAME(dcopy)(&nn, Vz, &incOne, Lz, &incOne);
    F77_NAME(dpotrf)(lower, &n, Lz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVy, &incOne);
    for(i = 0; i < n; i++){
      cholVy[i*n + i] += deltasq;
    }

    // find Cholesky of (Vz + deltasq*I)
    F77_NAME(dpotrf)(lower, &n, cholVy, &n, &info FCONE); if(info != 0){perror("c++ error: Vy dpotrf failed\n");}

    F77_NAME(dcopy)(&n, Y, &incOne, tmp_n1, &incOne);
    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);  // LyInv*y

    dtemp = pow(F77_NAME(dnrm2)(&n, tmp_n1, &incOne), 2);
    sse += dtemp;

    F77_NAME(dcopy)(&np, X, &incOne, tmp_np1, &incOne);
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, cholVy, &n, tmp_np1, &n FCONE FCONE FCONE FCONE); // LyInv*X
    F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np1, &n, tmp_n1, &incOne, &zero, tmp_p1, &incOne FCONE); // t(LyInv*y)*(LyInv*X)=Xt*VyInv*y

    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);   // Lbeta
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");}  // VbetaInv
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betaMu, &incOne, &zero, tmp_p2, &incOne FCONE);  // VbetaInv*muBeta

    dtemp = F77_CALL(ddot)(&p, betaMu, &incOne, tmp_p2, &incOne);  // t(muBeta)*VbetaInv*muBeta
    sse += dtemp;

    F77_NAME(daxpy)(&p, &one, tmp_p2, &incOne, tmp_p1, &incOne);  // Xt*Vyinv*y + VbetaInv*muBeta
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np1, &n, tmp_np1, &n, &zero, tmp_pp, &p FCONE FCONE);  // t(LyInv*X)*(LyInv*X)=Xt*VyInv*X
    F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, tmp_pp, &incOne);  // Xt*VyInv*X + VbetaInv

    F77_NAME(dcopy)(&p, tmp_p1, &incOne, tmp_p2, &incOne);
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}  // chol(Xt*VyInv*X + VbetaInv)
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p2, &incOne FCONE FCONE FCONE);  // chol(Xt*VyInv*X + VbetaInv)^{-1}(Xt*Vyinv*y + VbetaInv*muBeta)

    dtemp = pow(F77_NAME(dnrm2)(&p, tmp_p2, &incOne), 2);  // t(m)*M*m
    sse -= dtemp;

    // posterior parameters of sigmaSq
    sigmaSqIGaPost += sigmaSqIGa;
    sigmaSqIGaPost += 0.5 * n;

    sigmaSqIGbPost += sigmaSqIGb;
    sigmaSqIGbPost += 0.5 * sse;

    // set-up for sampling gamma = [beta, z]
    F77_NAME(dcopy)(&p, betaMu, &incOne, mu_vbeta, &incOne);                                   // mu_vbeta = muBeta
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, Lbeta, &p, mu_vbeta, &incOne FCONE FCONE FCONE);  // mu_vbeta = LbetaInv*muBeta

    // posterior samples of sigma-sq and beta
    SEXP samples_sigmaSq_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

    double sigmaSq = 0;
    double *v1 = (double *) R_alloc(p, sizeof(double)); zeros(v1, p);
    double *v2 = (double *) R_alloc(n, sizeof(double)); zeros(v2, n);
    double *out_p = (double *) R_alloc(p, sizeof(double)); zeros(out_p, p);
    double *out_n = (double *) R_alloc(n, sizeof(double)); zeros(out_n, n);

    GetRNGstate();

    for(s = 0; s < nSamples; s++){

      // sample sigmaSq from its marginal posterior
      dtemp = 1.0 / sigmaSqIGbPost;
      dtemp = rgamma(sigmaSqIGaPost, dtemp);
      sigmaSq = 1.0 / dtemp;
      REAL(samples_sigmaSq_r)[s] = sigmaSq;

      // sample fixed effects and spatial effects by projection
      dtemp = sqrt(sigmaSq);

      for(i = 0; i < n; i++){
        dtemp2 = deltaInv * Y[i];
        tmp_n1[i] = rnorm(dtemp2, dtemp);
        tmp_n1[i] = deltaInv * tmp_n1[i];
        tmp_n2[i] = rnorm(0.0, dtemp);
      }

      for(j = 0; j < p; j++){
        tmp_p1[j] = rnorm(mu_vbeta[j], dtemp);
      }

      F77_NAME(dtrsv)(lower, ytran, nUnit, &p, Lbeta, &p, tmp_p1, &incOne FCONE FCONE FCONE);
      F77_NAME(dtrsv)(lower, ytran, nUnit, &n, Lz, &n, tmp_n2, &incOne FCONE FCONE FCONE);

      F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n1, &incOne, &one, tmp_p1, &incOne FCONE);
      F77_NAME(daxpy)(&n, &one, tmp_n1, &incOne, tmp_n2, &incOne);

      F77_NAME(dcopy)(&p, tmp_p1, &incOne, v1, &incOne);
      F77_NAME(dcopy)(&n, tmp_n2, &incOne, v2, &incOne);

      inversionLM(X, n, p, deltasq, VbetaInv, Vz, cholVy, v1, v2,
                  tmp_n1, tmp_n2, tmp_p1, tmp_pp, tmp_np1,
                  out_p, out_n, 0);

      // inversionLM2(X, n, p, deltasq, VbetaInv, Vz, cholVy, v1, v2,
      //              out_p, out_n);

      F77_NAME(dcopy)(&p, &out_p[0], &incOne, &REAL(samples_beta_r)[s*p], &incOne);
      F77_NAME(dcopy)(&n, &out_n[0], &incOne, &REAL(samples_z_r)[s*n], &incOne);

    }

    PutRNGstate();

    // make return object for posterior samples of sigma-sq and beta
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

    // samples of sigma-sq
    SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

    // samples of z
    SET_VECTOR_ELT(result_r, 2, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z"));

    Rf_namesgets(result_r, resultName_r);

    // SEXP result_r = PROTECT(Rf_allocVector(REALSXP, nPp)); nProtect++;
    // SEXP result_r = PROTECT(Rf_allocMatrix(REALSXP, p, p)); nProtect++;
    // double *pointer_result_r = REAL(result_r);

    // for (i = 0; i < p; i++) {
    //   for (j = 0; j < n; j++) {
    //     REAL(result_r)[i*n + j] = tmp_np1[i*n + j];
    //   }
    // }

    // copyMatrixSEXP(tmp_pp, p, p, pointer_result_r);
    // copyVectorSEXP(v, nPp, pointer_result_r);

    UNPROTECT(nProtect);

    return result_r;

  }

}
