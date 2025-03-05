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
    char const *lside = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *nunit = "N";
    const double one = 1.0;
    const double negone = -1.0;
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

      if(loopd){
        if(loopd_method == exact_str){
          Rprintf("LOO-PD calculation method = %s\nNumber of Monte Carlo samples = %i.\n", loopd_method.c_str(), loopd_nMC);
        }
        if(loopd_method == cv_str){
          Rprintf("LOO-PD calculation method = %i-fold %s\nNumber of Monte Carlo samples = %i.\n", CV_K, loopd_method.c_str(), loopd_nMC);
        }
        if(loopd_method == psis_str){
          Rprintf("LOO-PD calculation method = %s\n", loopd_method.c_str());
        }
      }
      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/
    double dtemp1 = 0.0, dtemp2 = 0.0, dtemp3 = 0.0;

    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);                          // correlation matrix
    double *cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);                  // Cholesky of Vz
    double *cholVzPlusI = (double *) R_alloc(nn, sizeof(double)); zeros(cholVzPlusI, nn);        // allocate memory for n x n matrix
    double *cholSchur_n = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholSchur_n, nn);   // allocate memory for Schur complement
    double *cholSchur_p = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cholSchur_p, pp);   // allocate memory for Schur complement
    double *D1invX = (double *) R_chk_calloc(np, sizeof(double)); zeros(D1invX, np);             // allocate for preprocessing
    double *DinvB_pn = (double *) R_chk_calloc(np, sizeof(double)); zeros(DinvB_pn, np);         // allocate memory for p x n matrix
    double *DinvB_nn = (double *) R_chk_calloc(nn, sizeof(double)); zeros(DinvB_nn, nn);         // allocate memory for n x n matrix
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);              // allocate VbetaInv
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                    // Cholesky of Vbeta
    double *XtX = (double *) R_alloc(pp, sizeof(double)); zeros(XtX, pp);                        // Store XtX
    double *thetasp = (double *) R_alloc(2, sizeof(double));                                     // spatial process parameters

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

      SEXP loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

      // Exact leave-one-out predictive densities (LOO-PD) calculation
      if(loopd_method == exact_str){

        int n1 = n - 1;
        int n1n1 = n1 * n1;
        int n1p = n1 * p;

        // Set-up storage for leave-one-out data
        double *looY = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looY, n1);
        double *loo_nBinom = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_nBinom, n1);
        double *looX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looX, n1p);
        double *X_tilde = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_tilde, p);

        // Set-up storage for pre-processing
        double *looVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looVz, n1n1);
        double *looCholVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVz, n1n1);
        double *looCholVzPlusI = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVzPlusI, n1n1);
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

        // Set-up storage for spatial prediction
        double *looCz = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looCz, n1);
        double z_tilde = 0.0, z_tilde_var = 0.0, z_tilde_mu = 0.0, Mdist = 0.0;

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

          // Pre-processing for projGLM() on leave-one-out data
          cholRowDelUpdate(n, cholVz, loo_index, looCholVz, tmp_n11);                                            // Row-deletion CHOL update Vz
          cholRowDelUpdate(n, cholVzPlusI, loo_index, looCholVzPlusI, tmp_n11);                                  // Row-deletion CHOL update Vz+I

          // Spatial prediction variance term
          copyVecExcludingOne(&Vz[loo_index*n], looCz, n, loo_index);                                            // looCz = Vz[-i,i]
          F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholVz, &n1, looCz, &incOne FCONE FCONE FCONE);           // looCz = LzInv * Cz
          dtemp1 = pow(F77_NAME(dnrm2)(&n1, looCz, &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
          z_tilde_var = Vz[loo_index*n + loo_index] - dtemp1;                                                    // z_tilde_var = Vz_tilde - Czt*VzInv*Cz

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
            dtemp1 = 1.0 / dtemp3;
            dtemp2 = sqrt(dtemp1);
            for(j = 0; j < p; j++){
              loo_v_beta[j] = rnorm(0.0, dtemp2);                                                  // loo_v_beta ~ t_nu_beta(0, 1)
            }

            dtemp1 = 0.5 * nu_z;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp1 = 1.0 / dtemp3;
            dtemp2 = sqrt(dtemp1);
            for(loo_i = 0; loo_i < n1; loo_i++){
              loo_v_xi[loo_i] = rnorm(0.0, sigma_xi);                                              // loo_v_xi ~ N(0, 1)
              loo_v_z[loo_i] = rnorm(0.0, dtemp2);                                                 // loo_v_z ~ t_nu_z(0, 1)
            }

            // LOO projection step
            projGLM(looX, n1, p, loo_v_eta, loo_v_xi, loo_v_beta, loo_v_z, cholSchur_p1, cholSchur_n1, sigmaSq_xi, Lbeta,
                    looCholVz, looVz, looCholVzPlusI, D1invlooX, DinvB_pn1, DinvB_n1n1, tmp_n11, loo_tmp_p);

            // predict z at the loo_index location
            F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholVz, &n1, loo_v_z, &incOne FCONE FCONE FCONE);    // loo_v_z = LzInv * v_z
            z_tilde_mu = F77_CALL(ddot)(&n1, looCz, &incOne, loo_v_z, &incOne);                               // z_tilde_mu = Czt*VzInv*v_z
            Mdist = pow(F77_NAME(dnrm2)(&n1, loo_v_z, &incOne), 2);                                           // Mdist = v_zt*VzInv*v_z

            // sample z_tilde
            dtemp1 = 0.5 * (nu_z + n1);
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp1 = 1.0 / dtemp3;
            dtemp2 = dtemp1 * (Mdist + nu_z) / (nu_z + n1);
            dtemp3 = sqrt(dtemp2);
            z_tilde = rnorm(0.0, dtemp3);
            z_tilde = z_tilde * sqrt(z_tilde_var);
            z_tilde = z_tilde + z_tilde_mu;

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

        PutRNGstate();

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

      // K-fold cross-validation for LOO-PD calculation
      if(loopd_method == cv_str){

        int *startsCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(startsCV, CV_K);
        int *endsCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(endsCV, CV_K);
        int *sizesCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(sizesCV, CV_K);

        mkCVpartition(n, CV_K, startsCV, endsCV, sizesCV);

        int nk = 0;         // nk = size of k-th partition
        int nknk = 0;
        // int nknk = 0;       // nknk = nk x nk
        int nnk = 0;        // nnk = (n - nk); size of k-th partition deleted data
        // int nnknk = 0;
        // int nnknnk = 0;     // nnknnk = (n - nk) x (n - nk)
        // int nkp = 0;        // nkp = (n - nk) x p
        // int nnkp = 0;       // nnkp = (n - nk) x p

        int nkmin = findMin(sizesCV, CV_K);
        int nkmax = findMax(sizesCV, CV_K);
        int nknkmax = nkmax * nkmax;
        int nnkmax = n - nkmin;
        int nnknnkmax = nnkmax * nnkmax;
        int nnkmaxnkmax = nnkmax * nkmax;
        int nkmaxp = nkmax * p;
        int nnkmaxp = nnkmax * p;

        // Set-up storage for cross-validation data
        double *cvY = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cvY, nnkmax);                 // Store block-deleted Y
        double *cv_nBinom = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_nBinom, nnkmax);     // Store block-deleted nBinom
        double *cvX = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(cvX, nnkmaxp);               // Store block-deleted X

        // Set-up storage for pre-processing
        double *cvVz = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvVz, nnknnkmax);                         // Store block-deleted Vz
        double *cvCholVz = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholVz, nnknnkmax);                 // Store block-deleted Cholesky update of Vz
        double *cvCholVzPlusI = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholVzPlusI, nnknnkmax);       // Store block-deleted Chlesky update of Vz+I
        double *cvXtX = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cvXtX, pp);                                     // Store XtX
        double *DinvB_pnnkmax = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(DinvB_pnnkmax, nnkmaxp);           // allocate memory for p x max(n-nk) matrix
        double *DinvB_nnknnkmax = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(DinvB_nnknnkmax, nnknnkmax);   // allocate memory for max(n-nk) x max(n-nk) matrix
        double *cholSchur_p2 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cholSchur_p2, pp);                       // allocate memory for Schur complement
        double *cholSchur_nnkmax = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cholSchur_nnkmax, nnknnkmax); // allocate memory for Schur complement
        double *D1invcvX = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(D1invcvX, nnkmaxp);                     // allocate for preprocessing
        double *tmp_nnknnkmax = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(tmp_nnknnkmax, nnknnkmax);       // allocate for max(n-nk) x max(n-nk) matrix
        double *tmp_nnkmaxp = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(tmp_nnkmaxp, nnkmaxp);               // allocate for max(n-nk) x p matrix
        double *tmp_nnkmax = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(tmp_nnkmax, nnkmax);                   // allocate for max(n-nk) x 1 vector

        // Set-up storage for sampling for block-deleted model
        double *cv_v_eta = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_v_eta, nnkmax);
        double *cv_v_xi = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_v_xi, nnkmax);
        double *cv_v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(cv_v_beta, p);
        double *cv_v_z = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_v_z, nnkmax);
        double *cv_tmp_p = (double *) R_chk_calloc(p, sizeof(double)); zeros(cv_tmp_p, p);                       // temporary p x 1 vector

        // Set-up storage for held-out
        double *X_tilde = (double *) R_chk_calloc(nkmaxp, sizeof(double)); zeros(X_tilde, nkmaxp);         // Store held-out X
        double *Y_tilde = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(Y_tilde, nkmax);           // Store held-out Y
        double *nBinom_tilde = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(nBinom_tilde, nkmax); // Store held-out Y

        // Set-up storage for spatial prediction
        double *LzInvCz_cv = (double *) R_chk_calloc(nnkmaxnkmax, sizeof(double)); zeros(LzInvCz_cv, nnkmaxnkmax);  // cross-covariance matrix max(n-nk)xmax(nk)
        double *z_tilde_cov = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(z_tilde_cov, nknkmax);              // held-out covariance matrix max(nk)xmax(nk)
        double *tmp_nknkmax = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(tmp_nknkmax, nknkmax);        // allocate for max(nk) x max(nk) matrix
        double *z_tilde_mu = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(z_tilde_mu, nkmax);
        double *z_tilde = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(z_tilde, nkmax);
        double *tmp_nkmax = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(tmp_nkmax, nkmax);
        double PCMdist = 0.0;

        int cv_index = 0;
        int start_index = 0;
        int end_index = 0;
        int cv_i = 0;
        int sMC_CV = 0;
        int loopd_nMC_nkmax = loopd_nMC * nkmax;
        double *loopd_val_MC_CV = (double *) R_chk_calloc(loopd_nMC_nkmax, sizeof(double)); zeros(loopd_val_MC_CV, loopd_nMC_nkmax);

        GetRNGstate();

        for(cv_index = 0; cv_index < CV_K; cv_index++){

          nk = sizesCV[cv_index];
          nknk = nk * nk;
          nnk = n - nk;
          // nnkp = nnk * p;
          start_index = startsCV[cv_index];
          end_index = endsCV[cv_index];

          // Block-deleted data
          copyVecExcludingBlock(Y, cvY, n, start_index, end_index);                                                 // Block-deleted Y
          copyVecExcludingBlock(nBinom, cv_nBinom, n, start_index, end_index);                                      // Block-deleted nBinom
          copyMatrixDelRowBlock(X, n, p, cvX, start_index, end_index);                                              // Block-deleted X
          copyMatrixDelRowColBlock(Vz, n, n, cvVz, start_index, end_index, start_index, end_index);                 // Block-deleted Vz

          // Held-out data
          copyMatrixRowBlock(X, n, p, X_tilde, start_index, end_index);                                             // Held-out X = X_tilde
          copyVecBlock(Y, Y_tilde, n, start_index, end_index);                                                      // Held-out Y = Y_tilde
          copyVecBlock(nBinom, nBinom_tilde, n, start_index, end_index);                                            // Held-out nBinom = nBinom_tilde

          // Block-deleted Cholesky updates
          cholBlockDelUpdate(n, cholVz, start_index, end_index, cvCholVz, tmp_nnknnkmax, tmp_nnkmax);
          cholBlockDelUpdate(n, cholVzPlusI, start_index, end_index, cvCholVzPlusI, tmp_nnknnkmax, tmp_nnkmax);

          // Spatial process prediction
          copyMatrixRowColBlock(Vz, n, n, z_tilde_cov, start_index, end_index, start_index, end_index);                                 // z_tilde_cov = Vz[ids, ids]
          copyMatrixColDelRowBlock(Vz, n, n, LzInvCz_cv, start_index, end_index, start_index, end_index);                               // LzInvCz_cv = Vz[-ids, ids]
          F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, cvCholVz, &nnk, LzInvCz_cv, &nnk FCONE FCONE FCONE FCONE);       // LzInvCz_cv = inv(Lz)*Cz
          F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, LzInvCz_cv, &nnk, LzInvCz_cv, &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE); // tmp_nknkmax = t(Cz)*inv(Vz)*Cz
          F77_NAME(daxpy)(&nknk, &negone, tmp_nknkmax, &incOne, z_tilde_cov, &incOne);                                                  // z_tilde_cov = VzTilde - t(Cz)*inv(Vz)*Cz
          F77_NAME(dpotrf)(lower, &nk, z_tilde_cov, &nk, &info FCONE); if(info != 0){perror("c++ error: z_tilde_schur dpotrf failed\n");}
          mkLT(z_tilde_cov, nk);

          // Pre-processing for projGLM() on block-deleted data
          F77_NAME(dgemm)(ytran, ntran, &p, &p, &nnk, &one, cvX, &nnk, cvX, &nnk, &zero, cvXtX, &p FCONE FCONE);    // XtX = t(X)*X
          cholSchurGLM(cvX, nnk, p, sigmaSq_xi, cvXtX, VbetaInv, cvVz, cvCholVzPlusI, tmp_nnknnkmax, tmp_nnkmaxp,
                       DinvB_pnnkmax, DinvB_nnknnkmax, cholSchur_p2, cholSchur_nnkmax, D1invcvX);

          // Fit on block-deleted data and obtain LOO-PD by Monte Carlo average
          for(sMC_CV = 0; sMC_CV < loopd_nMC; sMC_CV++){

            if(family == family_poisson){
              for(cv_i = 0; cv_i < nnk; cv_i++){
                dtemp1 = cvY[cv_i] + epsilon;
                dtemp2 = 1.0;
                dtemp3 = rgamma(dtemp1, dtemp2);
                cv_v_eta[cv_i] = log(dtemp3);
              }
            }

            if(family == family_binomial){
              for(cv_i = 0; cv_i < nnk; cv_i++){
                dtemp1 = cvY[cv_i] + epsilon;
                dtemp2 = cv_nBinom[cv_i];
                dtemp2 += 2.0 * epsilon;
                dtemp2 -= dtemp1;
                dtemp3 = rbeta(dtemp1, dtemp2);
                cv_v_eta[cv_i] = logit(dtemp3);
              }
            }

            if(family == family_binary){
              for(cv_i = 0; cv_i < nnk; cv_i++){
                dtemp1 = cvY[cv_i] + epsilon;
                dtemp2 = cv_nBinom[cv_i];
                dtemp2 += 2.0 * epsilon;
                dtemp2 -= dtemp1;
                dtemp3 = rbeta(dtemp1, dtemp2);
                cv_v_eta[cv_i] = logit(dtemp3);
              }
            }

            dtemp1 = 0.5 * nu_beta;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp3 = sqrt(dtemp3);
            for(j = 0; j < p; j++){
              cv_v_beta[j] = rnorm(0.0, dtemp3);                                                  // loo_v_beta ~ N(0, 1)
            }

            dtemp1 = 0.5 * nu_z;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp3 = sqrt(dtemp3);
            for(cv_i = 0; cv_i < nnk; cv_i++){
              cv_v_xi[cv_i] = rnorm(0.0, sigma_xi);                                              // loo_v_xi ~ N(0, 1)
              cv_v_z[cv_i] = rnorm(0.0, dtemp3);                                                 // loo_v_z ~ N(0, 1)
            }

            // LOO projection step
            projGLM(cvX, nnk, p, cv_v_eta, cv_v_xi, cv_v_beta, cv_v_z, cholSchur_p2, cholSchur_nnkmax, sigmaSq_xi, Lbeta,
                    cvCholVz, cvVz, cvCholVzPlusI, D1invcvX, DinvB_pnnkmax, DinvB_nnknnkmax, tmp_nnkmax, cv_tmp_p);

            // Prediction of spatial process at held-out locations
            F77_NAME(dtrsv)(lower, ntran, nunit, &nnk, cvCholVz, &nnk, cv_v_z, &incOne FCONE FCONE FCONE);                // cv_v_z = LzInv * v_z
            F77_NAME(dgemv)(ytran, &nnk, &nk, &one, LzInvCz_cv, &nnk, cv_v_z, &incOne, &zero, z_tilde_mu, &incOne FCONE); // z_tilde_mu = t(Cz)*inv(Vz)*v_z
            PCMdist = pow(F77_NAME(dnrm2)(&nnk, cv_v_z, &incOne), 2);                                                     // Mahalanobis distance t(v_z)&inv(Vz)*v_z

            // sample z_tilde
            dtemp1 = 0.5 * (nu_z + nnk);
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp2 = 1.0 / dtemp3;
            dtemp1 = (PCMdist + nu_z) / (nu_z + nnk);
            dtemp3 = dtemp1 * dtemp2;
            dtemp1 = sqrt(dtemp3);
            for(cv_i = 0; cv_i < nk; cv_i++){
              z_tilde[cv_i] = rnorm(0.0, dtemp1);
            }
            F77_NAME(dgemv)(ntran, &nk, &nk, &one, z_tilde_cov, &nk, z_tilde, &incOne, &zero, tmp_nkmax, &incOne FCONE);
            F77_NAME(daxpy)(&nk, &one, z_tilde_mu, &incOne, tmp_nkmax, &incOne);
            F77_NAME(dcopy)(&nk, tmp_nkmax, &incOne, z_tilde, &incOne);

            // Find canonical parameter (X_tilde*v_beta + z_tilde)
            F77_NAME(dgemv)(ntran, &nk, &p, &one, X_tilde, &nk, cv_v_beta, &incOne, &zero, tmp_nkmax, &incOne FCONE);
            F77_NAME(daxpy)(&nk, &one, z_tilde, &incOne, tmp_nkmax, &incOne);

            // Find CV-LOO-PD
            if(family == family_poisson){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = exp(tmp_nkmax[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dpois(Y_tilde[cv_i], dtemp1, 1);
              }
            }

            if(family == family_binomial){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = inverse_logit(tmp_nkmax[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dbinom(Y_tilde[cv_i], nBinom_tilde[cv_i], dtemp1, 1);
              }
            }

            if(family == family_binary){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = inverse_logit(tmp_nkmax[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dbinom(Y_tilde[cv_i], 1.0, dtemp1, 1);
              }
            }

          }

          for(cv_i = 0; cv_i < nk; cv_i++){
            REAL(loopd_out_r)[start_index + cv_i] = logMeanExp(&loopd_val_MC_CV[cv_i*loopd_nMC], loopd_nMC);
          }
        }

        PutRNGstate();

        R_chk_free(startsCV);
        R_chk_free(endsCV);
        R_chk_free(sizesCV);
        R_chk_free(cvY);
        R_chk_free(cvX);
        R_chk_free(cv_nBinom);
        R_chk_free(cvVz);
        R_chk_free(cvCholVz);
        R_chk_free(cvCholVzPlusI);
        R_chk_free(cvXtX);
        R_chk_free(DinvB_pnnkmax);
        R_chk_free(DinvB_nnknnkmax);
        R_chk_free(cholSchur_p2);
        R_chk_free(cholSchur_nnkmax);
        R_chk_free(D1invcvX);
        R_chk_free(tmp_nnkmax);
        R_chk_free(tmp_nnkmaxp);
        R_chk_free(tmp_nnknnkmax);
        R_chk_free(cv_v_eta);
        R_chk_free(cv_v_xi);
        R_chk_free(cv_v_beta);
        R_chk_free(cv_v_z);
        R_chk_free(cv_tmp_p);
        R_chk_free(X_tilde);
        R_chk_free(Y_tilde);
        R_chk_free(nBinom_tilde);
        R_chk_free(LzInvCz_cv);
        R_chk_free(z_tilde_cov);
        R_chk_free(tmp_nknkmax);
        R_chk_free(z_tilde_mu);
        R_chk_free(z_tilde);
        R_chk_free(tmp_nkmax);
        R_chk_free(loopd_val_MC_CV);

      }

      // Pareto-smoothed Importance Sampling for LOO-PD calculation
      if(loopd_method == psis_str){

        int loo_index = 0, s = 0;
        double theta_i = 0.0, z_s = 0.0;
        double dtemp_psis;

        double *X_i = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_i, p);
        double *beta_s = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta_s, p);

        double *dens_i = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(dens_i, nSamples);
        double *rawIR = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(rawIR, nSamples);
        double *sortedIR = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(sortedIR, nSamples);
        double *stableIR = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(stableIR, nSamples);
        int *orderIR = (int *) R_chk_calloc(nSamples, sizeof(int)); zeros(orderIR, nSamples);

        // find M = floor(min(0.2*S, 3*sqrt(S)))
        double val1 = 0.0, val2 = 0.0, min_val = 0.0;
        int M = 0;
        val1 = 0.2 * nSamples;
        val2 = 3 * sqrt(nSamples);
        min_val = fmin2(val1, val2);
        M = (int)floor(min_val);

        double *tmp_M1 = (double *) R_chk_calloc(M, sizeof(double)); zeros(tmp_M1, M);
        double *tmp_M2 = (double *) R_chk_calloc(M, sizeof(double)); zeros(tmp_M2, M);
        double *tmp_M3 = (double *) R_chk_calloc(M, sizeof(double)); zeros(tmp_M3, M);
        double *ksigma = (double *) R_chk_calloc(2, sizeof(double)); zeros(ksigma, 2);

        double *pointer_beta = REAL(samples_beta_r);
        double *pointer_z = REAL(samples_z_r);

        for(loo_index = 0; loo_index < n; loo_index++){

          copyMatrixRowToVec(X, n, p, X_i, loo_index);                          // X_i = X[i,1:p]

          for(s = 0; s < nSamples; s++){

            copyMatrixColToVec(pointer_beta, p, nSamples, beta_s, s);           // beta_s = beta[s], s-th sample
            z_s = pointer_z[n*s + loo_index];                                   // z_s = z_i[s], s-th sample of i-th spatial effect
            theta_i = F77_CALL(ddot)(&p, X_i, &incOne, beta_s, &incOne);        // theta_i = X_i * beta_s
            theta_i += z_s;                                                     // theta_i = X_i*beta_s + zi_s

            if(family == family_poisson){
              dtemp_psis = exp(theta_i);
              dens_i[s] = dpois(Y[loo_index], dtemp_psis, 1);
              rawIR[s] = - dens_i[s];
            }

            if(family == family_binomial){
              dtemp_psis = inverse_logit(theta_i);
              dens_i[s] = dbinom(Y[loo_index], nBinom[loo_index], dtemp_psis, 1);
              rawIR[s] = - dens_i[s];
            }

            if(family == family_binary){
              dtemp_psis = inverse_logit(theta_i);
              dens_i[s] = dbinom(Y[loo_index], 1.0, dtemp_psis, 1);
              rawIR[s] = - dens_i[s];
            }

          }

          ParetoSmoothedIR(rawIR, M, nSamples, sortedIR, orderIR, stableIR, ksigma, tmp_M1, tmp_M2, tmp_M3);

          REAL(loopd_out_r)[loo_index] = logWeightedSumExp(dens_i, rawIR, nSamples);

        }

        R_chk_free(X_i);
        R_chk_free(beta_s);
        R_chk_free(dens_i);
        R_chk_free(rawIR);
        R_chk_free(sortedIR);
        R_chk_free(stableIR);
        R_chk_free(orderIR);
        R_chk_free(tmp_M1);
        R_chk_free(tmp_M2);
        R_chk_free(tmp_M3);
        R_chk_free(ksigma);

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