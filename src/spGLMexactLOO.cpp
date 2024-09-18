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

  SEXP spGLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP nBinom_r,
                     SEXP coordsD_r, SEXP corfn_r, SEXP betaV_r, SEXP nu_beta_r,
                     SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP phi_r, SEXP nu_r,
                     SEXP epsilon_r, SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                     SEXP CV_K_r, SEXP loopd_nMC_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *nUnit = "N";
    const double one = 1.0;
    const double zero = 0.0;
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

    double *coordsD = REAL(coordsD_r);

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

    // Leave-one-out predictive density details
    int loopd = INTEGER(loopd_r)[0];
    std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));
    int CV_K = INTEGER(CV_K_r)[0];
    int loopd_nMC = INTEGER(loopd_nMC_r)[0];

    const char *exact_str = "exact";
    const char *cv_str = "cv";
    const char *psis_str = "psis";

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
      Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
      Rprintf("\n");

      Rprintf("\tsigmaSq.beta ~ IG(nu.beta/2, nu.beta/2), nu.beta = %.2f\n", nu_beta);
      Rprintf("\tsigmaSq.z ~ IG(nu.z/2, nu.z/2),\tnu.z = %.2f\n", nu_z);
      Rprintf("\tsigmaSq.xi = %.2f\n", sigmaSq_xi);
      Rprintf("\tBoundary adjustment parameter = %.2f\n\n", epsilon);

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.3f, and, nu = %.3f\n\n", phi, nu);
      }else{
        Rprintf("\tphi = %.3f\n\n", phi);
      }

      Rprintf("Number of posterior samples = %i.\n", nSamples);

      if(loopd){

        if(loopd_method == exact_str){
          Rprintf("Finding leave-one-out predictive densities (LOO-PD) using\n");
          Rprintf("method = %s, and, number of Monte Carlo samples = %i.\n", loopd_method.c_str(), loopd_nMC);
        }

        if(loopd_method == cv_str){
          Rprintf("Finding leave-one-out predictive densities (LOO-PD) using\n");
          Rprintf("method = %i-fold %s, and, number of Monte Carlo samples = %i.\n", CV_K, loopd_method.c_str(), loopd_nMC);
        }

        if(loopd_method == psis_str){
          Rprintf("Finding leave-one-out predictive densities (LOO-PD) using\n");
          Rprintf("method = %s, (Pareto-smoothed Importance Sampling)", loopd_method.c_str());
        }

      }

      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/
    double dtemp1, dtemp2, dtemp3;

    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);                       // correlation matrix
    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);               // Cholesky of Vz
    double *cholVzPlusI = (double *) R_alloc(nn, sizeof(double)); zeros(cholVzPlusI, nn);     // allocate memory for n x n matrix
    double *cholSchur_n = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholSchur_n, nn);     // allocate memory for Schur complement
    double *cholSchur_p = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cholSchur_p, pp);     // allocate memory for Schur complement
    double *D1invX = (double *) R_chk_calloc(np, sizeof(double)); zeros(D1invX, np);               // allocate for preprocessing
    double *DinvB_pn = (double *) R_chk_calloc(np, sizeof(double)); zeros(DinvB_pn, np);           // allocate memory for p x n matrix
    double *DinvB_nn = (double *) R_chk_calloc(nn, sizeof(double)); zeros(DinvB_nn, nn);           // allocate memory for n x n matrix
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);           // allocate VbetaInv
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                 // Cholesky of Vbeta
    double *XtX = (double *) R_alloc(pp, sizeof(double)); zeros(XtX, pp);                     // Store XtX
    double *thetasp = (double *) R_alloc(2, sizeof(double));                                  // spatial process parameters

    //construct covariance matrix (full)
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull(coordsD, n, thetasp, corfn, Vz);

    // Find Cholesky of Vz
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}

    // construct unit spherical perturbation of Vz; (Vz+I)
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVzPlusI, &incOne);
    for(i = 0; i < n; i++){
      cholVzPlusI[i*n + i] += 1.0;
    }

    // find Cholesky factor of unit spherical perturbation of Vz
    F77_NAME(dpotrf)(lower, &n, cholVzPlusI, &n, &info FCONE); if(info != 0){perror("c++ error: VzPlusI dpotrf failed\n");}

    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                     // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);                                                     // Lbeta = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");} // VbetaInv = chol2inv(Vbeta)

    // Find XtX
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, X, &n, &zero, XtX, &p FCONE FCONE);                   // XtX = t(X)*X

    // Get the Schur complement of top left nxn submatrix of (HtH)
    double *tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np, np);       // temporary allocate memory for n x p matrix
    double *tmp_nn = (double *) R_chk_calloc(nn, sizeof(double)); zeros(tmp_nn, nn);       // temporary allocate memory for n x n matrix

    cholSchurGLM(X, n, p, sigmaSq_xi, XtX, VbetaInv, Vz, cholVzPlusI, tmp_nn, tmp_np,
                 DinvB_pn, DinvB_nn, cholSchur_p, cholSchur_n, D1invX);

    R_chk_free(tmp_nn);
    R_chk_free(tmp_np);

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

    double *v_eta = (double *) R_chk_calloc(n, sizeof(double)); zeros(v_eta, n);
    double *v_xi = (double *) R_chk_calloc(n, sizeof(double)); zeros(v_xi, n);
    double *v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(v_beta, p);
    double *v_z = (double *) R_chk_calloc(n, sizeof(double)); zeros(v_z, n);

    double *tmp_n = (double *) R_chk_calloc(n, sizeof(double)); zeros(tmp_n, n);           // allocate memory for n x 1 vector
    double *tmp_p = (double *) R_chk_calloc(p, sizeof(double)); zeros(tmp_p, p);           // allocate memory for p x 1 vector

    GetRNGstate();

    for(s = 0; s < nSamples; s++){

      if(family == family_poisson){
        for(i = 0; i < n; i++){
          dtemp1 = Y[i] + epsilon;
          dtemp2 = 1.0;
          dtemp3 = rgamma(dtemp1, dtemp2);
          v_eta[i] = log(dtemp3);
        }
      }

      if(family == family_binomial){
        for(i = 0; i < n; i++){
          dtemp1 = Y[i] + epsilon;
          dtemp2 = nBinom[i];
          dtemp2 += 2.0 * epsilon;
          dtemp2 -= dtemp1;
          dtemp3 = rbeta(dtemp1, dtemp2);
          v_eta[i] = logit(dtemp3);
        }
      }

      if(family == family_binary){
        for(i = 0; i < n; i++){
          dtemp1 = Y[i] + epsilon;
          dtemp2 = nBinom[i];
          dtemp2 += 2.0 * epsilon;
          dtemp2 -= dtemp1;
          dtemp3 = rbeta(dtemp1, dtemp2);
          v_eta[i] = logit(dtemp3);
        }
      }

      dtemp1 = 0.5 * nu_beta;
      dtemp2 = 1.0 / dtemp1;
      dtemp3 = rgamma(dtemp1, dtemp2);
      dtemp3 = 1.0 / dtemp3;
      dtemp3 = sqrt(dtemp3);
      for(j = 0; j < p; j++){
        v_beta[j] = rnorm(0.0, dtemp3);                                                  // v_beta ~ N(0, 1)
      }

      dtemp1 = 0.5 * nu_z;
      dtemp2 = 1.0 / dtemp1;
      dtemp3 = rgamma(dtemp1, dtemp2);
      dtemp3 = 1.0 / dtemp3;
      dtemp3 = sqrt(dtemp3);
      for(i = 0; i < n; i++){
        v_xi[i] = rnorm(0.0, sigma_xi);                                                  // v_xi ~ N(0, 1)
        v_z[i] = rnorm(0.0, dtemp3);                                                     // v_z ~ N(0, 1)
      }

      // projection step
      projGLM(X, n, p, v_eta, v_xi, v_beta, v_z, cholSchur_p, cholSchur_n, sigmaSq_xi, Lbeta,
              cholVz, Vz, cholVzPlusI, D1invX, DinvB_pn, DinvB_nn, tmp_n, tmp_p);

      // copy samples into SEXP return object
      F77_NAME(dcopy)(&p, &v_beta[0], &incOne, &REAL(samples_beta_r)[s*p], &incOne);
      F77_NAME(dcopy)(&n, &v_z[0], &incOne, &REAL(samples_z_r)[s*n], &incOne);
      F77_NAME(dcopy)(&n, &v_xi[0], &incOne, &REAL(samples_xi_r)[s*n], &incOne);

    }

    PutRNGstate();

    R_chk_free(tmp_n);
    R_chk_free(tmp_p);
    R_chk_free(v_eta);
    R_chk_free(v_xi);
    R_chk_free(v_beta);
    R_chk_free(v_z);
    R_chk_free(cholSchur_n);
    R_chk_free(cholSchur_p);
    R_chk_free(D1invX);
    R_chk_free(DinvB_pn);
    R_chk_free(DinvB_nn);

    // make return object
    SEXP result_r, resultName_r;

    if(loopd){

      int n1 = n - 1;
      int n1n1 = n1 * n1;
      int n1p = n1 * p;

      SEXP loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

      // Exact leave-one-out predictive densities (LOO-PD) calculation
      if(loopd_method == exact_str){

        // Set-up storage for pre-processing
        double *looY = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looY, n1);
        double *loo_nBinom = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_nBinom, n1);
        double *looX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looX, n1p);
        double *X_tilde = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_tilde, p);
        double *looVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looVz, n1n1);
        double *looCholVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVz, n1n1);
        double *looCholVzPlusI = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVzPlusI, n1n1);
        double *looCz = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looCz, n1);
        double *looXtX = (double *) R_chk_calloc(pp, sizeof(double)); zeros(looXtX, pp);                           // Store XtX
        double *DinvB_pn1 = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(DinvB_pn1, n1p);                   // allocate memory for p x n matrix
        double *DinvB_n1n1 = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(DinvB_n1n1, n1n1);               // allocate memory for n x n matrix
        double *cholSchur_n1 = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(cholSchur_n1, n1n1);           // allocate memory for Schur complement
        double *cholSchur_p1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cholSchur_p1, pp);               // allocate memory for Schur complement
        double *D1invlooX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(D1invlooX, n1p);                   // allocate for preprocessing
        double *tmp_n11 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n11, n1);

        // Get the Schur complement of top left n1xn1 submatrix of (HtH)
        double *tmp_n1p = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(tmp_n1p, n1p);                       // temporary n1 x p matrix
        double *tmp_n1n1 = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(tmp_n1n1, n1n1);                   // temporary n1 x n1 matrix

        // Set-up storage for sampling for leave-one-out model fit
        double *loo_v_eta = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_v_eta, n1);
        double *loo_v_xi = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_v_xi, n1);
        double *loo_v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(loo_v_beta, p);
        double *loo_v_z = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_v_z, n1);
        double *loo_tmp_p = (double *) R_chk_calloc(p, sizeof(double)); zeros(loo_tmp_p, p);                       // temporary p x 1 vector
        double z_tilde, z_tilde_var, z_tilde_mu;

        int loo_index = 0;
        int loo_i = 0;
        int sMC = 0;
        double *loopd_val_MC = (double *) R_chk_calloc(loopd_nMC, sizeof(double)); zeros(loopd_val_MC, loopd_nMC);

        GetRNGstate();

        for(loo_index = 0; loo_index < n; loo_index++){

          // Prepare leave-one-out data
          copyVecExcludingOne(Y, looY, n, loo_index);                                                            // Leave-one-out Y
          copyVecExcludingOne(nBinom, loo_nBinom, n, loo_index);                                                 // Leave-one-out nBinom
          copyMatrixDelRow(X, n, p, looX, loo_index);                                                            // Row-deleted X
          copyMatrixRowToVec(X, n, p, X_tilde, loo_index);                                                       // Copy left out X into Xtilde
          copyMatrixDelRowCol(Vz, n, n, looVz, loo_index, loo_index);                                            // Row-column deleted Vz
          copyVecExcludingOne(&Vz[loo_index*n], looCz, n, loo_index);                                            // looCz = Vz[-i,i]

          // Pre-processing for projGLM() on leave-one-out data
          cholRowDelUpdate(n, cholVz, loo_index, looCholVz, tmp_n11);                                            // Row-deletion CHOL update Vz
          cholRowDelUpdate(n, cholVzPlusI, loo_index, looCholVzPlusI, tmp_n11);                                  // Row-deletion CHOL update Vy
          F77_NAME(dgemm)(ytran, ntran, &p, &p, &n1, &one, looX, &n1, looX, &n1, &zero, looXtX, &p FCONE FCONE); // XtX = t(X)*X
          cholSchurGLM(looX, n1, p, sigmaSq_xi, looXtX, VbetaInv, looVz, looCholVzPlusI, tmp_n1n1, tmp_n1p,
                       DinvB_pn1, DinvB_n1n1, cholSchur_p1, cholSchur_n1, D1invlooX);

          for(sMC = 0; sMC < loopd_nMC; sMC++){

            if(family == family_poisson){
              for(loo_i = 0; loo_i < n1; loo_i++){
                dtemp1 = looY[loo_i] + epsilon;
                dtemp2 = 1.0;
                dtemp3 = rgamma(dtemp1, dtemp2);
                loo_v_eta[loo_i] = log(dtemp3);
              }
            }

            if(family == family_binomial){
              for(loo_i = 0; loo_i < n1; loo_i++){
                dtemp1 = looY[loo_i] + epsilon;
                dtemp2 = loo_nBinom[loo_i];
                dtemp2 += 2.0 * epsilon;
                dtemp2 -= dtemp1;
                dtemp3 = rbeta(dtemp1, dtemp2);
                loo_v_eta[loo_i] = logit(dtemp3);
              }
            }

            if(family == family_binary){
              for(loo_i = 0; loo_i < n1; loo_i++){
                dtemp1 = looY[loo_i] + epsilon;
                dtemp2 = loo_nBinom[loo_i];
                dtemp2 += 2.0 * epsilon;
                dtemp2 -= dtemp1;
                dtemp3 = rbeta(dtemp1, dtemp2);
                loo_v_eta[loo_i] = logit(dtemp3);
              }
            }

            dtemp1 = 0.5 * nu_beta;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp3 = sqrt(dtemp3);
            for(j = 0; j < p; j++){
              loo_v_beta[j] = rnorm(0.0, dtemp3);                                                  // loo_v_beta ~ N(0, 1)
            }

            dtemp1 = 0.5 * nu_z;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp3 = sqrt(dtemp3);
            for(loo_i = 0; loo_i < n1; loo_i++){
              loo_v_xi[loo_i] = rnorm(0.0, sigma_xi);                                              // loo_v_xi ~ N(0, 1)
              loo_v_z[loo_i] = rnorm(0.0, dtemp3);                                                 // loo_v_z ~ N(0, 1)
            }

            // LOO projection step
            projGLM(looX, n1, p, loo_v_eta, loo_v_xi, loo_v_beta, loo_v_z, cholSchur_p1, cholSchur_n1, sigmaSq_xi, Lbeta,
                    looCholVz, looVz, looCholVzPlusI, D1invlooX, DinvB_pn1, DinvB_n1n1, tmp_n11, loo_tmp_p);

            // predict z at the loo_index location
            F77_NAME(dcopy)(&n1, looCz, &incOne, tmp_n11, &incOne);
            F77_NAME(dtrsv)(lower, ntran, nUnit, &n1, looCholVz, &n1, tmp_n11, &incOne FCONE FCONE FCONE);    // tmp_n11 = LzInv * Cz
            F77_NAME(dtrsv)(lower, ntran, nUnit, &n1, looCholVz, &n1, loo_v_z, &incOne FCONE FCONE FCONE);    // loo_v_z = LzInv * v_z
            z_tilde_var = Vz[loo_index*n + loo_index];
            dtemp1 = pow(F77_NAME(dnrm2)(&n1, tmp_n11, &incOne), 2);                                          // dtemp1 = Czt*VzInv*Cz
            z_tilde_var -= dtemp1;                                                                            // z_tilde_var = VzTilde - Czt*VzInv*Cz
            dtemp1 = pow(F77_NAME(dnrm2)(&n1, loo_v_z, &incOne), 2);                                          // dtemp1 = v_zt*VzInv*v_z
            dtemp1 += nu_z;                                                                                   // dtemp1 = nu_z + v_zt*VzInv*v_z
            dtemp2 = dtemp1 / (nu_z + n1);                                                                    // dtemp2 = (nu_z+v_zt*VzInv*v_z)/(nu_z+n1)
            dtemp3 = dtemp2 * z_tilde_var;
            z_tilde_var = dtemp3;                                                                             // z_tilde_var = dtemp2*(VzTilde - Czt*VzInv*Cz)
            z_tilde_mu = F77_CALL(ddot)(&p, tmp_n11, &incOne, loo_v_z, &incOne);                              // z_tilde_mu = Czt*VzInv*v_z

            // sample z_tilde
            dtemp1 = 0.5 * (nu_z + n1);
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp1 = dtemp3 * z_tilde_var;
            dtemp2 = sqrt(dtemp1);
            z_tilde = rnorm(z_tilde_mu, dtemp2);

            dtemp1 = F77_CALL(ddot)(&p, X_tilde, &incOne, loo_v_beta, &incOne);
            dtemp2 = dtemp1 + z_tilde;                                                                        // dtemp2 = X_tilde*beta + z_tilde

            // Find predictive densities from canonical parameter dtemp2 = (X*beta + z)
            if(family == family_poisson){
              dtemp3 = exp(dtemp2);
              loopd_val_MC[sMC] = dpois(Y[loo_index], dtemp3, 1);
            }

            if(family == family_binomial){
              dtemp3 = inverse_logit(dtemp2);
              loopd_val_MC[sMC] = dbinom(Y[loo_index], nBinom[loo_index], dtemp3, 1);
            }

            if(family == family_binary){
              dtemp3 = inverse_logit(dtemp2);
              loopd_val_MC[sMC] = dbinom(Y[loo_index], 1.0, dtemp3, 1);
            }

          }

          REAL(loopd_out_r)[loo_index] = logMeanExp(loopd_val_MC, loopd_nMC);

        }

        R_chk_free(looY);
        R_chk_free(loo_nBinom);
        R_chk_free(looX);
        R_chk_free(X_tilde);
        R_chk_free(looVz);
        R_chk_free(looCholVz);
        R_chk_free(looCholVzPlusI);
        R_chk_free(looXtX);
        R_chk_free(DinvB_pn1);
        R_chk_free(DinvB_n1n1);
        R_chk_free(cholSchur_n1);
        R_chk_free(cholSchur_p1);
        R_chk_free(D1invlooX);
        R_chk_free(tmp_n11);
        R_chk_free(tmp_n1p);
        R_chk_free(tmp_n1n1);
        R_chk_free(loo_v_eta);
        R_chk_free(loo_v_xi);
        R_chk_free(loo_v_beta);
        R_chk_free(loo_v_z);
        R_chk_free(loo_tmp_p);

      }

      PutRNGstate();

      // K-fold cross-validation for LOO-PD calculation
      if(loopd_method == cv_str){

        int loo_index = 0;

        for(loo_index = 0; loo_index < n; loo_index++){
          REAL(loopd_out_r)[loo_index] = 0.0;
        }

      }

      // Pareto-smoothed Importance Sampling for LOO-PD calculation
      if(loopd_method == psis_str){

        int loo_index = 0;

        for(loo_index = 0; loo_index < n; loo_index++){
          REAL(loopd_out_r)[loo_index] = 0.0;
        }

      }

      // make return object for posterior samples and leave-one-out predictive densities
      int nResultListObjs = 4;

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

      // loo-pd
      // leave-one-out predictive densities
      SET_VECTOR_ELT(result_r, 3, loopd_out_r);
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("loopd"));

      Rf_namesgets(result_r, resultName_r);

    }else{

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

      // samples of xi
      SET_VECTOR_ELT(result_r, 2, samples_xi_r);
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("xi"));

      Rf_namesgets(result_r, resultName_r);

    }



    UNPROTECT(nProtect);

    return result_r;

  } // end spGLMexact
}