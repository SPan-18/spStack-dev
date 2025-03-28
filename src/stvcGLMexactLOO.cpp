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

  SEXP stvcGLMexactLOO(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                       SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                       SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                       SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                       SEXP CV_K_r, SEXP loopd_nMC_r,  SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, k, s, info, nProtect = 0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *nunit = "N";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    double *Y = REAL(Y_r);
    double *nBinom = REAL(nBinom_r);
    double *X = REAL(X_r);
    double *X_tilde = REAL(X_tilde_r);
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int np = n * p;
    int r = INTEGER(r_r)[0];
    int rr = r * r;
    int nr = n * r;
    int nrp = nr * p;
    int nnr = nn * r;
    int nrnr = nr * nr;

    std::string family = CHAR(STRING_ELT(family_r, 0));

    double *coords_sp = REAL(sp_coords_r);
    double *coords_tm = REAL(time_coords_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // priors
    double *betaMu = (double *) R_alloc(p, sizeof(double)); zeros(betaMu, p);
    double *betaV = (double *) R_alloc(pp, sizeof(double)); zeros(betaV, pp);
    F77_NAME(dcopy)(&pp, REAL(betaV_r), &incOne, betaV, &incOne);

    double nu_beta = REAL(nu_beta_r)[0];
    double nu_z = REAL(nu_z_r)[0];
    double sigmaSq_xi = REAL(sigmaSq_xi_r)[0];
    double sigma_xi = sqrt(sigmaSq_xi);

    double *iwScale = (double *) R_alloc(rr, sizeof(double)); zeros(iwScale, rr);
    F77_NAME(dcopy)(&rr, REAL(iwScale_r), &incOne, iwScale, &incOne);

    // spatial-temporal process parameters: create spatial-temporal covariance matrices
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    double *Vz = NULL;
    double *R = NULL;

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

        }else if(processType == "multivariate2"){

          phi_s_vec[0] = REAL(phi_s_r)[0];
          phi_t_vec[0] = REAL(phi_t_r)[0];
          thetaspt[0] = phi_s_vec[0];
          thetaspt[1] = phi_t_vec[0];

          R = (double *) R_alloc(nn, sizeof(double)); zeros(R, nn);
          sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, R);
          Vz = (double *) R_alloc(nrnr, sizeof(double)); zeros(Vz, nrnr);
          kronecker(r, n, iwScale, R, Vz);

        }
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
      Rprintf("\tMODEL DESCRIPTION\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Family = %s.\n\n", family.c_str());
      Rprintf("Number of fixed effects = %i.\n", p);
      Rprintf("Number of varying coefficients = %i.\n\n", r);

      Rprintf("Priors:\n");

      Rprintf("\tbeta: Gaussian\n");
      Rprintf("\tmu:"); printVec(betaMu, p);
      Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
      Rprintf("\n");

      Rprintf("\tsigmaSq.beta ~ IG(nu.beta/2, nu.beta/2)\n");
      Rprintf("\tnu.beta = %.2f, nu.z = %.2f.\n", nu_beta, nu_z);
      Rprintf("\tSpatial-temporal process model: %s.\n", processType.c_str());
      if(processType == "multivariate"){
        Rprintf("\tSigma: Inverse-Wishart\n");
        Rprintf("\tdf: %.2f\n", nu_z);
        Rprintf("\tScale:\n"); printMtrx(iwScale, r, r);
      }else{
        Rprintf("\tsigmaSq.z.j ~ IG(nu.z/2, nu.z/2), j = 1,...,%i.\n", r);
      }
      Rprintf("\tsigmaSq.xi = %.2f.\n", sigmaSq_xi);
      Rprintf("\tBoundary adjustment parameter = %.2f.\n\n", epsilon);

      Rprintf("Spatial-temporal correlation function: %s.\n", corfn.c_str());

      Rprintf("Process type: %s.\n", processType.c_str());

      if(processType == "independent.shared" || processType == "multivariate"){
        Rprintf("All %i spatial-temporal processes share common parameters:\n", r);
        if(corfn == "gneiting-decay"){
            Rprintf("\tphi_s = %.2f, and, phi_t = %.2f.\n\n", phi_s_vec[0], phi_t_vec[0]);
        }
      }else{
        Rprintf("Parameters for the %i spatial-temporal process(es):\n", r);
        if(corfn == "gneiting-decay"){
            Rprintf("\tphi_s ="); printVec(phi_s_vec, r);
            Rprintf("\tphi_t ="); printVec(phi_t_vec, r);
        }
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

    double *cholVz = NULL;               // define NULL pointer for chol(Vz)
    double *cholR = NULL;
    double *chol_iwScale = NULL;

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
      chol_iwScale = (double *) R_alloc(rr, sizeof(double)); zeros(chol_iwScale, rr);
      F77_NAME(dcopy)(&rr, iwScale, &incOne, chol_iwScale, &incOne);
      F77_NAME(dpotrf)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){perror("c++ error: iwScale dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){perror("c++ error: iwScale dpotri failed\n");} // chol_iwScale = chol2inv(iwScale)
      F77_NAME(dpotrf)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){perror("c++ error: inv(iwScale) dpotrf failed\n");}
      mkLT(chol_iwScale, r);

    }else if(processType == "multivariate2"){

        // // Inefficient cholesky!
        // cholVz = (double *) R_alloc(nrnr, sizeof(double)); zeros(cholVz, nrnr);          // nrxnr matrix chol(Vz)
        // F77_NAME(dcopy)(&nrnr, Vz, &incOne, cholVz, &incOne);
        // F77_NAME(dpotrf)(lower, &nr, cholVz, &nr, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
        // mkLT(cholVz, nr);

        // Efficient Cholesky for kronecker product: chol(kron(A, B)) = kron(chol(A), chol(B))
        cholVz = (double *) R_alloc(nrnr, sizeof(double)); zeros(cholVz, nrnr);          // nrxnr matrix chol(Vz)
        cholR = (double *) R_alloc(nn, sizeof(double)); zeros(cholR, nn);
        chol_iwScale = (double *) R_alloc(rr, sizeof(double)); zeros(chol_iwScale, rr);
        F77_NAME(dcopy)(&nn, R, &incOne, cholR, &incOne);
        F77_NAME(dcopy)(&rr, iwScale, &incOne, chol_iwScale, &incOne);
        F77_NAME(dpotrf)(lower, &n, cholR, &n, &info FCONE); if(info != 0){perror("c++ error: R dpotrf failed\n");}
        F77_NAME(dpotrf)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){perror("c++ error: iwScale dpotrf failed\n");}
        chol_kron(r, n, chol_iwScale, cholR, cholVz);
        mkLT(cholVz, nr);

    }

    // Allocations for XtX, XTildetX, and VbetaInv
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);           // allocate VbetaInv
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                 // Cholesky of Vbeta
    double *XtX = (double *) R_alloc(pp, sizeof(double)); zeros(XtX, pp);                     // Store XtX
    double *XTildetX = (double *) R_alloc(nrp, sizeof(double)); zeros(XTildetX, nrp);         // Store XTildetX

    // Find VbetaInv
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                           // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: VBeta dpotrf failed\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);                                                           // Lbeta = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");}       // VbetaInv = chol2inv(Vbeta)

    // Find XtX
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, X, &n, &zero, XtX, &p FCONE FCONE);                   // XtX = t(X)*X

    // Find t(X_tilde)*X
    lmulm_XTilde_VC(ytran, n, r, p, X_tilde, X, XTildetX);

    // Allocations for I + Xtilde*Vz*t(Xtilde)
    double *XTildeVzXTildet = (double *) R_alloc(nn, sizeof(double)); zeros(XTildeVzXTildet, nn);
    double *cholIplusXTildeVzXTildet = (double *) R_alloc(nn, sizeof(double)); zeros(cholIplusXTildeVzXTildet, nn);
    double *VzXTildet = (double *) R_chk_calloc(nnr, sizeof(double)); zeros(VzXTildet, nnr);

    rmul_Vz_XTildeT(n, r, X_tilde, Vz, VzXTildet, processType);                                                  // Vz*t(X_tilde)
    lmulm_XTilde_VC(ntran, n, r, n, X_tilde, VzXTildet, XTildeVzXTildet);                                        // X_tilde*Vz*t(X_tilde)
    R_chk_free(VzXTildet);

    // Find Cholesky of capacitance matric: I + XTilde*Vz*t(XTilde)
    F77_NAME(dcopy)(&nn, XTildeVzXTildet, &incOne, cholIplusXTildeVzXTildet, &incOne);
    for(i = 0; i < n; i++){
        cholIplusXTildeVzXTildet[i*n + i] += 1.0;
    }
    F77_NAME(dpotrf)(lower, &n, cholIplusXTildeVzXTildet, &n, &info FCONE);
    if(info != 0){perror("c++ error: capacitance matrix dpotrf failed\n");}
    mkLT(cholIplusXTildeVzXTildet, n);

    // Allocations for priming step (pre-processing)
    double *tmp_nnr = (double *) R_chk_calloc(nnr, sizeof(double)); zeros(tmp_nnr, nnr);
    double *D1Inv = (double *) R_chk_calloc(nrnr, sizeof(double)); zeros(D1Inv, nrnr);
    double *D1InvB1 = (double *) R_chk_calloc(nrp, sizeof(double)); zeros(D1InvB1, nrp);
    double *cholschurA1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cholschurA1, pp);
    double *DInvB_pn = (double *) R_chk_calloc(np, sizeof(double)); zeros(DInvB_pn, np);
    double *DInvB_nrn = (double *) R_chk_calloc(nnr, sizeof(double)); zeros(DInvB_nrn, nnr);
    double *cholschurA = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholschurA, nn);
    double *tmp_rr = (double *) R_alloc(rr, sizeof(double)); zeros(tmp_rr, rr);
    double *samp_Sigma = (double *) R_alloc(rr, sizeof(double)); zeros(samp_Sigma, rr);

    // Evaluate priming step
    primingGLMvc(n, p, r, X, X_tilde, XtX, XTildetX, VbetaInv, Vz, processType, cholIplusXTildeVzXTildet,
                 sigmaSq_xi, tmp_nnr, D1Inv, D1InvB1, cholschurA1, DInvB_pn, DInvB_nrn, cholschurA);

    R_chk_free(tmp_nnr);

    /*****************************************
     Set-up posterior sampling
     *****************************************/
    // posterior samples of sigma-sq and beta
    SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, nr, nSamples)); nProtect++;
    SEXP samples_xi_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

    const char *family_poisson = "poisson";
    const char *family_binary = "binary";
    const char *family_binomial = "binomial";

    double *v_eta = (double *) R_chk_calloc(n, sizeof(double)); zeros(v_eta, n);
    double *v_xi = (double *) R_chk_calloc(n, sizeof(double)); zeros(v_xi, n);
    double *v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(v_beta, p);
    double *v_z = (double *) R_chk_calloc(nr, sizeof(double)); zeros(v_z, nr);
    double *tmp_nr = (double *) R_chk_calloc(nr, sizeof(double)); zeros(tmp_nr, nr);

    double dtemp1 = 0.0, dtemp2 = 0.0, dtemp3 = 0.0;

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

      for(i = 0; i < n; i++){
        v_xi[i] = rnorm(0.0, sigma_xi);                                                  // v_xi ~ N(0, sigmaSq_xi)
      }

      dtemp1 = 0.5 * nu_beta;
      dtemp2 = 1.0 / dtemp1;
      dtemp3 = rgamma(dtemp1, dtemp2);
      dtemp3 = 1.0 / dtemp3;
      dtemp3 = sqrt(dtemp3);
      for(j = 0; j < p; j++){
        v_beta[j] = rnorm(0.0, dtemp3);                                                  // v_beta ~ t
      }

      if(processType == "independent.shared"){
        dtemp1 = 0.5 * nu_z;
        dtemp2 = 1.0 / dtemp1;
        dtemp3 = rgamma(dtemp1, dtemp2);
        dtemp3 = 1.0 / dtemp3;
        dtemp3 = sqrt(dtemp3);
        for(k = 0; k < r; k++){
          for(i = 0; i < n; i++){
            v_z[k*n + i] = rnorm(0.0, dtemp3);
          }
        }
      }else if(processType == "independent"){
        for(k = 0; k < r; k++){
          dtemp1 = 0.5 * nu_z;
          dtemp2 = 1.0 / dtemp1;
          dtemp3 = rgamma(dtemp1, dtemp2);
          dtemp3 = 1.0 / dtemp3;
          dtemp3 = sqrt(dtemp3);
          for(i = 0; i < n; i++){
            v_z[k*n + i] = rnorm(0.0, dtemp3);
          }
        }
      }else if(processType == "multivariate"){

        for(k = 0; k < r; k++){
          for(i = 0; i < n; i++){
            tmp_nr[k*n + i] = rnorm(0.0, 1.0);
          }
        }
        rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr);
        F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){perror("c++ error: samp_Sigma dpotrf failed\n");}
        F77_NAME(dgemm)(ntran, ytran, &n, &r, &r, &one, tmp_nr, &n, samp_Sigma, &r, &zero, v_z, &n FCONE FCONE);

      }

      // projection step
      projGLMvc(n, p, r, X, X_tilde, sigmaSq_xi, Lbeta, cholVz, processType,
                v_eta, v_xi, v_beta, v_z, D1Inv, D1InvB1, cholschurA1,
                DInvB_pn, DInvB_nrn, cholschurA, tmp_nr);

      // copy samples into SEXP return object
      F77_NAME(dcopy)(&p, &v_beta[0], &incOne, &REAL(samples_beta_r)[s*p], &incOne);
      F77_NAME(dcopy)(&nr, &v_z[0], &incOne, &REAL(samples_z_r)[s*nr], &incOne);
      F77_NAME(dcopy)(&n, &v_xi[0], &incOne, &REAL(samples_xi_r)[s*n], &incOne);

    }

    PutRNGstate();

    R_chk_free(v_eta);
    R_chk_free(v_xi);
    R_chk_free(v_beta);
    R_chk_free(v_z);
    R_chk_free(tmp_nr);

    R_chk_free(D1Inv);
    R_chk_free(D1InvB1);
    R_chk_free(cholschurA1);
    R_chk_free(DInvB_pn);
    R_chk_free(DInvB_nrn);
    R_chk_free(cholschurA);

    // make return object
    SEXP result_r, resultName_r;

    if(loopd){

      if(verbose){
        Rprintf("Evaluating leave-one-out predictive densities.\n");
      }

      SEXP loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

      // Exact leave-one-out predictive densities (LOO-PD) calculation
      if(loopd_method == exact_str){

        int n1 = n - 1;
        int n1n1 = n1 * n1;
        int n1p = n1 * p;
        int n1r = n1 * r;
        int n1rp = n1r * p;
        int n1n1r = n1n1 * r;
        int n1rn1r = n1r * n1r;

        // Set-up storage for leave-one-out data
        double *looY = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looY, n1);
        double *loo_nBinom = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_nBinom, n1);
        double *looX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looX, n1p);
        double *looX_tilde = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(looX_tilde, n1r);
        double *X_pred = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_pred, p);
        double *X_tilde_pred = (double *) R_chk_calloc(r, sizeof(double)); zeros(X_tilde_pred, r);

        // Set-up storage for pre-processing for leave-one-out data
        double *looVz = NULL;
        double *looCholVz = NULL;
        double *looCholR = NULL;
        double *looCz = NULL;

        if(corfn == "gneiting-decay"){

          if(processType == "independent.shared" || processType == "multivariate"){

            looVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looVz, n1n1);
            looCholVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVz, n1n1);
            looCz = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looCz, n1);

          }else if(processType == "independent"){

            looVz = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(looVz, n1n1r);
            looCholVz = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(looCholVz, n1n1r);
            looCz = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(looCz, n1r);

          }else if(processType == "multivariate2"){

            looVz = (double *) R_chk_calloc(n1rn1r, sizeof(double)); zeros(looVz, n1rn1r);
            looCholVz = (double *) R_chk_calloc(n1rn1r, sizeof(double)); zeros(looCholVz, n1rn1r);
            looCholR = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholR, n1n1);
            looCz = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looCz, n1);

          }

        }

        double *looXtX = (double *) R_chk_calloc(pp, sizeof(double)); zeros(looXtX, pp);
        double *looXTildetX = (double *) R_chk_calloc(n1rp, sizeof(double)); zeros(looXTildetX, n1rp);
        double *looCholIplusXTildeVzXTildet = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholIplusXTildeVzXTildet, n1n1);

        // set-up pre-processing memory allocations for priming on leave-one-out data
        double *looD1Inv = (double *) R_chk_calloc(n1rn1r, sizeof(double)); zeros(looD1Inv, n1rn1r);
        double *looD1InvB1 = (double *) R_chk_calloc(n1rp, sizeof(double)); zeros(looD1InvB1, n1rp);
        double *looCholschurA1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(looCholschurA1, pp);
        double *looDInvB_pn = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looDInvB_pn, n1p);
        double *looDInvB_nrn = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(looDInvB_nrn, n1n1r);
        double *looCholschurA = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholschurA, n1n1);
        double *tmp_n1n1r = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(tmp_n1n1r, n1n1r);
        double *tmp_n11 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n11, n1);
        double *tmp_n1r = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(tmp_n1r, n1r);

        // Set-up storage for sampling for leave-one-out model fit
        double *loo_v_eta = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_v_eta, n1);
        double *loo_v_xi = (double *) R_chk_calloc(n1, sizeof(double)); zeros(loo_v_xi, n1);
        double *loo_v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(loo_v_beta, p);
        double *loo_v_z = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(loo_v_z, n1r);

        int loo_index = 0;
        int loo_i = 0;
        int sMC = 0;
        double *loopd_val_MC = (double *) R_chk_calloc(loopd_nMC, sizeof(double)); zeros(loopd_val_MC, loopd_nMC);
        double *z_tilde_var = (double *) R_chk_calloc(r, sizeof(double)); zeros(z_tilde_var, r);
        double *z_tilde_mu = (double *) R_chk_calloc(r, sizeof(double)); zeros(z_tilde_mu, r);
        double *z_tilde = (double *) R_chk_calloc(r, sizeof(double)); zeros(z_tilde, r);
        double *Mdist_r = (double *) R_chk_calloc(r, sizeof(double)); zeros(Mdist_r, r);
        double *Mdist_rr = (double *) R_chk_calloc(rr, sizeof(double)); zeros(Mdist_rr, rr);
        double *tmp_r = (double *) R_chk_calloc(r, sizeof(double)); zeros(tmp_r, r);

        GetRNGstate();

        for(loo_index = 0; loo_index < n; loo_index++){

          // Prepare leave-one-out data
          copyVecExcludingOne(Y, looY, n, loo_index);                           // Leave-one-out Y
          copyVecExcludingOne(nBinom, loo_nBinom, n, loo_index);                // Leave-one-out nBinom
          copyMatrixDelRow(X, n, p, looX, loo_index);                           // Row-deleted X
          copyMatrixDelRow(X_tilde, n, r, looX_tilde, loo_index);               // Row-deleted X_tilde
          copyMatrixRowToVec(X, n, p, X_pred, loo_index);                       // Copy left out X into X_pred
          copyMatrixRowToVec(X_tilde, n, r, X_tilde_pred, loo_index);           // Copy left out X into X_pred

          // Leave-one-out XtX, substract Xi*t(Xi) from XtX, instead of multiplying again
          F77_NAME(dgemm)(ntran, ytran, &p, &p, &incOne, &one, X_pred, &p, X_pred, &p, &zero, looXtX, &p FCONE FCONE);
          F77_NAME(dscal)(&pp, &negOne, looXtX, &incOne);
          F77_NAME(daxpy)(&pp, &one, XtX, &incOne, looXtX, &incOne);

          // Leave-one-out XTildetX, delete i-th row from every n-th block of original XTildetX, instead of multiplying again
          copyMatrixDelRow_vc(XTildetX, nr, p, looXTildetX, loo_index, n);

          // Constructing leave-one-out Vz for each spatial-temporal process model, and also Schur complement for prediction
          if(processType == "independent.shared" || processType == "multivariate"){

            copyMatrixDelRowCol(Vz, n, n, looVz, loo_index, loo_index);
            cholRowDelUpdate(n, cholVz, loo_index, looCholVz, tmp_n11);

            copyVecExcludingOne(&Vz[loo_index*n], looCz, n, loo_index);                                            // looCz = Vz[-i,i]
            F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholVz, &n1, looCz, &incOne FCONE FCONE FCONE);           // looCz = LzInv * Cz
            dtemp1 = pow(F77_NAME(dnrm2)(&n1, looCz, &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
            z_tilde_var[0] = Vz[loo_index*n + loo_index] - dtemp1;                                                 // z_tilde_var = Vz_tilde - Czt*VzInv*Cz

          }else if(processType == "independent"){

            for(k = 0; k < r; k++){
              copyMatrixDelRowCol(&Vz[nn*k], n, n, &looVz[n1n1*k], loo_index, loo_index);
              cholRowDelUpdate(n, &cholVz[nn*k], loo_index, &looCholVz[n1n1*k], tmp_n11);

              copyVecExcludingOne(&Vz[nn*k + loo_index*n], &looCz[n1*k], n, loo_index);                                     // looCz = Vz[-i,i]
              F77_NAME(dtrsv)(lower, ntran, nunit, &n1, &looCholVz[n1n1*k], &n1, &looCz[n1*k], &incOne FCONE FCONE FCONE);  // looCz = LzInv * Cz
              dtemp1 = pow(F77_NAME(dnrm2)(&n1, &looCz[n1*k], &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
              z_tilde_var[k] = Vz[nn*k + loo_index*n + loo_index] - dtemp1;

            }

          }else if(processType == "multivariate2"){

            copyMatrixDelRowCol_vc(Vz, nr, nr, looVz, loo_index, loo_index, n);   // Row-column deleted Vz
            cholRowDelUpdate(n, cholR, loo_index, looCholR, tmp_n11);             // Row-column deleted R
            chol_kron(r, n1, chol_iwScale, looCholR, looCholVz);                  // Find kron(chol(Psi), chol(looR))
            mkLT(looCholVz, n1r);

            copyVecExcludingOne(&R[loo_index*n], looCz, n, loo_index);                                             // looCz = R[-i,i]
            F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholR, &n1, looCz, &incOne FCONE FCONE FCONE);            // looCz = LzInv * Cz
            dtemp1 = pow(F77_NAME(dnrm2)(&n1, looCz, &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
            z_tilde_var[0] = R[loo_index*n + loo_index] - dtemp1;                                                  // z_tilde_var = Vz_tilde - Czt*VzInv*Cz

          }

          // Constructing leave-one-out Cholesky factor for I + XTilde*Vz*t(XTilde)
          // It can be shown that it is equivalent to updating Cholesky factor
          // after removing i-th row and i-th column from original chol(I + XTilde*Vz*t(XTilde))
          cholRowDelUpdate(n, cholIplusXTildeVzXTildet, loo_index, looCholIplusXTildeVzXTildet, tmp_n11);

          primingGLMvc(n1, p, r, looX, looX_tilde, looXtX, looXTildetX, VbetaInv, looVz, processType, looCholIplusXTildeVzXTildet,
                       sigmaSq_xi, tmp_n1n1r, looD1Inv, looD1InvB1, looCholschurA1, looDInvB_pn, looDInvB_nrn, looCholschurA);

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

            for(loo_i = 0; loo_i < n1; loo_i++){
              loo_v_xi[loo_i] = rnorm(0.0, sigma_xi);
            }

            dtemp1 = 0.5 * nu_beta;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp1 = 1.0 / dtemp3;
            dtemp2 = sqrt(dtemp1);
            for(j = 0; j < p; j++){
              loo_v_beta[j] = rnorm(0.0, dtemp2);                                                  // loo_v_beta ~ t_nu_beta(0, 1)
            }

            if(processType == "independent.shared"){
              dtemp1 = 0.5 * nu_z;
              dtemp2 = 1.0 / dtemp1;
              dtemp3 = rgamma(dtemp1, dtemp2);
              dtemp3 = 1.0 / dtemp3;
              dtemp3 = sqrt(dtemp3);
              for(k = 0; k < r; k++){
                for(loo_i = 0; loo_i < n1; loo_i++){
                  loo_v_z[k*n1 + loo_i] = rnorm(0.0, dtemp3);
                }
              }
            }else if(processType == "independent"){
              for(k = 0; k < r; k++){
                dtemp1 = 0.5 * nu_z;
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp3 = 1.0 / dtemp3;
                dtemp3 = sqrt(dtemp3);
                for(loo_i = 0; loo_i < n1; loo_i++){
                  loo_v_z[k*n1 + loo_i] = rnorm(0.0, dtemp3);
                }
              }
            }else if(processType == "multivariate"){
              for(k = 0; k < r; k++){
                for(loo_i = 0; loo_i < n1; loo_i++){
                  tmp_n1r[k*n1 + loo_i] = rnorm(0.0, 1.0);
                }
              }
              rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr);
              F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){perror("c++ error: samp_Sigma dpotrf failed\n");}
              mkLT(samp_Sigma, r);
              F77_NAME(dgemm)(ntran, ytran, &n1, &r, &r, &one, tmp_n1r, &n1, samp_Sigma, &r, &zero, loo_v_z, &n1 FCONE FCONE);

            }

            // projection step
            projGLMvc(n1, p, r, looX, looX_tilde, sigmaSq_xi, Lbeta, looCholVz, processType,
                      loo_v_eta, loo_v_xi, loo_v_beta, loo_v_z, looD1Inv, looD1InvB1, looCholschurA1,
                      looDInvB_pn, looDInvB_nrn, looCholschurA, tmp_n1r);

            // Prediction at held-out point for each spatial-temporal process model
            if(processType == "independent.shared"){

              for(k = 0; k < r; k++){
                F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholVz, &n1, &loo_v_z[n1*k], &incOne FCONE FCONE FCONE);    // loo_v_z = LzInv * v_z
                z_tilde_mu[k] = F77_CALL(ddot)(&n1, looCz, &incOne, &loo_v_z[n1*k], &incOne);                            // z_tilde_mu = Czt*VzInv*v_z
                Mdist_r[k] = pow(F77_NAME(dnrm2)(&n1, &loo_v_z[n1*k], &incOne), 2);                                      // Mdist = v_zt*VzInv*v_z

                // sample z_tilde
                dtemp1 = 0.5 * (nu_z + n1);
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp1 = 1.0 / dtemp3;
                dtemp2 = dtemp1 * (Mdist_r[k] + nu_z) / (nu_z + n1);
                dtemp3 = sqrt(dtemp2);
                z_tilde[k] = rnorm(0.0, dtemp3);
                z_tilde[k] = z_tilde[k] * sqrt(z_tilde_var[0]);
                z_tilde[k] = z_tilde[k] + z_tilde_mu[k];
              }

            }else if(processType == "independent"){

                for(k = 0; k < r; k++){
                  F77_NAME(dtrsv)(lower, ntran, nunit, &n1, &looCholVz[n1n1*k], &n1, &loo_v_z[n1*k], &incOne FCONE FCONE FCONE); // loo_v_z = LzInv * v_z
                  z_tilde_mu[k] = F77_CALL(ddot)(&n1, &looCz[n1*k], &incOne, &loo_v_z[n1*k], &incOne);                           // z_tilde_mu = Czt*VzInv*v_z
                  Mdist_r[k] = pow(F77_NAME(dnrm2)(&n1, &loo_v_z[n1*k], &incOne), 2);                                            // Mdist = v_zt*VzInv*v_z

                  // sample z_tilde
                  dtemp1 = 0.5 * (nu_z + n1);
                  dtemp2 = 1.0 / dtemp1;
                  dtemp3 = rgamma(dtemp1, dtemp2);
                  dtemp1 = 1.0 / dtemp3;
                  dtemp2 = dtemp1 * (Mdist_r[k] + nu_z) / (nu_z + n1);
                  dtemp3 = sqrt(dtemp2);
                  z_tilde[k] = rnorm(0.0, dtemp3);
                  z_tilde[k] = z_tilde[k] * sqrt(z_tilde_var[k]);
                  z_tilde[k] = z_tilde[k] + z_tilde_mu[k];
              }

            }else if(processType == "multivariate"){

              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n1, &r, &one, looCholVz, &n1, loo_v_z, &n1 FCONE FCONE FCONE FCONE);          // loo_v_z = cholinv(Vz)*Z
              F77_NAME(dgemm)(ytran, ntran, &incOne, &r, &n1, &one, looCz, &n1, loo_v_z, &n1, &zero, z_tilde_mu, &incOne FCONE FCONE);   // z_tilde_mu = t(C)*inv(R)*Z
              F77_NAME(dgemm)(ytran, ntran, &r, &r, &n1, &one, loo_v_z, &n1, loo_v_z, &n1, &zero, Mdist_rr, &r FCONE FCONE);             // Mdist = t(Z)*inv(R)*Z
              F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, Mdist_rr, &incOne);                                                           // Mdist = iwScale + t(Z)*inv(R)*Z
              F77_NAME(dpotrf)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              F77_NAME(dpotri)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotri failed\n");}
              F77_NAME(dpotrf)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              mkLT(Mdist_rr, r);
              rInvWishart(r, nu_z + n1 + 2*r, Mdist_rr, samp_Sigma, tmp_rr);
              F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){perror("c++ error: samp_Sigma dpotrf failed\n");}
              mkLT(samp_Sigma, r);

              dtemp1 = sqrt(z_tilde_var[0]);
              for(k = 0; k < r; k++){
                tmp_r[k] = rnorm(0.0, dtemp1);
              }
              F77_NAME(dgemv)(ntran, &r, &r, &one, samp_Sigma, &r, tmp_r, &incOne, &zero, z_tilde, &incOne FCONE);
              F77_NAME(daxpy)(&r, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

            }else if(processType == "multivariate2"){

              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n1, &r, &one, looCholR, &n1, loo_v_z, &n1 FCONE FCONE FCONE FCONE);          // loo_v_z = cholinv(R)*Z
              F77_NAME(dgemm)(ytran, ntran, &incOne, &r, &n1, &one, looCz, &n1, loo_v_z, &n1, &zero, z_tilde_mu, &incOne FCONE FCONE);  // z_tilde_mu = t(C)*inv(R)*Z
              F77_NAME(dgemm)(ytran, ntran, &r, &r, &n1, &one, loo_v_z, &n1, loo_v_z, &n1, &zero, Mdist_rr, &r FCONE FCONE);            // Mdist = t(Z)*inv(R)*Z
              F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, Mdist_rr, &incOne);                                                          // Mdist = iwScale + t(Z)*inv(R)*Z
              F77_NAME(dpotrf)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              mkLT(Mdist_rr, r);

              // sample z_tilde
              dtemp1 = 0.5 * (nu_z + n1);
              dtemp2 = 1.0 / dtemp1;
              dtemp3 = rgamma(dtemp1, dtemp2);
              dtemp1 = 1.0 / dtemp3;
              dtemp3 = sqrt(dtemp1 * z_tilde_var[0]);
              for(k = 0; k < r; k++){
                z_tilde[k] = rnorm(0.0, dtemp3);
              }
              // multiply chol(Mdist) * z_tilde, then add z_tilde_mu
              F77_NAME(dgemv)(ntran, &r, &r, &one, Mdist_rr, &r, z_tilde, &incOne, &zero, tmp_r, &incOne FCONE);
              F77_NAME(dcopy)(&r, tmp_r, &incOne, z_tilde, &incOne);
              F77_NAME(daxpy)(&r, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

            }

            dtemp1 = F77_CALL(ddot)(&p, X_pred, &incOne, loo_v_beta, &incOne);
            dtemp1 += F77_CALL(ddot)(&r, X_tilde_pred, &incOne, z_tilde, &incOne);

            // Find predictive densities from canonical parameter dtemp2 = (X*beta + z)
            if(family == family_poisson){
              dtemp2 = exp(dtemp1);
              loopd_val_MC[sMC] = dpois(Y[loo_index], dtemp2, 1);
            }

            if(family == family_binomial){
              dtemp2 = inverse_logit(dtemp1);
              loopd_val_MC[sMC] = dbinom(Y[loo_index], nBinom[loo_index], dtemp2, 1);
            }

            if(family == family_binary){
              dtemp2 = inverse_logit(dtemp1);
              loopd_val_MC[sMC] = dbinom(Y[loo_index], 1.0, dtemp2, 1);
            }

          }

          REAL(loopd_out_r)[loo_index] = logMeanExp(loopd_val_MC, loopd_nMC);

        }

        PutRNGstate();

        R_chk_free(looY);
        R_chk_free(loo_nBinom);
        R_chk_free(looX);
        R_chk_free(looX_tilde);
        R_chk_free(X_pred);
        R_chk_free(X_tilde_pred);
        R_chk_free(looVz);
        R_chk_free(looCholVz);
        R_chk_free(looCholR);
        R_chk_free(looCz);
        R_chk_free(looXtX);
        R_chk_free(looXTildetX);
        R_chk_free(looCholIplusXTildeVzXTildet);
        R_chk_free(looD1Inv);
        R_chk_free(looD1InvB1);
        R_chk_free(looCholschurA1);
        R_chk_free(looDInvB_pn);
        R_chk_free(looDInvB_nrn);
        R_chk_free(looCholschurA);
        R_chk_free(tmp_n1n1r);
        R_chk_free(tmp_n11);
        R_chk_free(tmp_n1r);
        R_chk_free(loo_v_eta);
        R_chk_free(loo_v_xi);
        R_chk_free(loo_v_beta);
        R_chk_free(loo_v_z);
        R_chk_free(loopd_val_MC);
        R_chk_free(z_tilde_var);
        R_chk_free(z_tilde_mu);
        R_chk_free(z_tilde);
        R_chk_free(Mdist_r);
        R_chk_free(Mdist_rr);
        R_chk_free(tmp_r);

      }

      // K-fold cross-validation for LOO-PD calculation
      if(loopd_method == cv_str){

        int *startsCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(startsCV, CV_K);
        int *endsCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(endsCV, CV_K);
        int *sizesCV = (int *) R_chk_calloc(CV_K, sizeof(int)); zeros(sizesCV, CV_K);

        mkCVpartition(n, CV_K, startsCV, endsCV, sizesCV);

        int nk = 0;         // nk = size of k-th partition
        int nknk = 0;
        int nnk = 0;
        int nnknnk = 0;
        int nnkr = 0;
        int nkr = 0;

        int nkmin = findMin(sizesCV, CV_K);
        int nkmax = findMax(sizesCV, CV_K);
        int nknkmax = nkmax * nkmax;
        int nnkmax = n - nkmin;
        int nnknnkmax = nnkmax * nnkmax;
        int nnknnkmaxr = nnknnkmax * r;
        int nnkmaxnkmax = nnkmax * nkmax;
        int nnkmaxnkmaxr = nnkmaxnkmax * r;
        int nknkmaxr = nknkmax * r;
        int nkmaxp = nkmax * p;
        int nkmaxr = nkmax * r;
        int nnkmaxp = nnkmax * p;
        int nnkmaxr = nnkmax*r;
        int nnkmaxrp = nnkmaxr * p;
        int nnkmaxrnnkmaxr = nnkmaxr * nnkmaxr;

        // Set-up storage for cross-validation data
        double *cvY = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cvY, nnkmax);                   // Store block-deleted Y
        double *cv_nBinom = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_nBinom, nnkmax);       // Store block-deleted nBinom
        double *cvX = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(cvX, nnkmaxp);                 // Store block-deleted X
        double *cvX_tilde = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(cvX_tilde, nnkmaxr);     // Store block-deleted X
        double *X_pred = (double *) R_chk_calloc(nkmaxp, sizeof(double)); zeros(X_pred, nkmaxp);             // Store held-out X
        double *X_tilde_pred = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(X_tilde_pred, nkmaxr); // Store held-out X
        double *Y_pred = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(Y_pred, nkmax);               // Store held-out Y
        double *nBinom_pred = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(nBinom_pred, nkmax);     // Store held-out X

        // Set-up storage for pre-processing for cross-validated data
        double *cvVz = NULL;
        double *cvCholVz = NULL;
        double *cvCholR = NULL;
        double *cvCz = NULL;
        double *z_tilde_cov = NULL;
        double *z_tilde_mu = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(z_tilde_mu, nkmaxr);
        double *z_tilde = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(z_tilde, nkmaxr);
        double *PCM_dist = (double *) R_chk_calloc(rr, sizeof(double)); zeros(PCM_dist, rr);

        if(corfn == "gneiting-decay"){

          if(processType == "independent.shared" || processType == "multivariate"){

            cvVz = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvVz, nnknnkmax);
            cvCholVz = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholVz, nnknnkmax);
            cvCz = (double *) R_chk_calloc(nnkmaxnkmax, sizeof(double)); zeros(cvCz, nnkmaxnkmax);
            z_tilde_cov = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(z_tilde_cov, nknkmax);

          }else if(processType == "independent"){

            cvVz = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(cvVz, nnknnkmaxr);
            cvCholVz = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(cvCholVz, nnknnkmaxr);
            cvCz = (double *) R_chk_calloc(nnkmaxnkmaxr, sizeof(double)); zeros(cvCz, nnkmaxnkmaxr);
            z_tilde_cov = (double *) R_chk_calloc(nknkmaxr, sizeof(double)); zeros(z_tilde_cov, nknkmaxr);

          }else if(processType == "multivariate2"){

            cvVz = (double *) R_chk_calloc(nnkmaxrnnkmaxr, sizeof(double)); zeros(cvVz, nnkmaxrnnkmaxr);
            cvCholVz = (double *) R_chk_calloc(nnkmaxrnnkmaxr, sizeof(double)); zeros(cvCholVz, nnkmaxrnnkmaxr);
            cvCholR = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholR, nnknnkmax);
            cvCz = (double *) R_chk_calloc(nnkmaxnkmax, sizeof(double)); zeros(cvCz, nnkmaxnkmax);
            z_tilde_cov = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(z_tilde_cov, nknkmax);

          }

        }

        double *cvXtX = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cvXtX, pp);
        double *cvXTildetX = (double *) R_chk_calloc(nnkmaxrp, sizeof(double)); zeros(cvXTildetX, nnkmaxrp);
        double *cvCholIplusXTildeVzXTildet = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholIplusXTildeVzXTildet, nnknnkmax);

        // set-up pre-processing memory allocations for priming on leave-one-out data
        double *cvD1Inv = (double *) R_chk_calloc(nnkmaxrnnkmaxr, sizeof(double)); zeros(cvD1Inv, nnkmaxrnnkmaxr);
        double *cvD1InvB1 = (double *) R_chk_calloc(nnkmaxrp, sizeof(double)); zeros(cvD1InvB1, nnkmaxrp);
        double *cvCholschurA1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cvCholschurA1, pp);
        double *cvDInvB_pn = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(cvDInvB_pn, nnkmaxp);
        double *cvDInvB_nrn = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(cvDInvB_nrn, nnknnkmaxr);
        double *cvCholschurA = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholschurA, nnknnkmax);
        double *tmp_n1n1r = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(tmp_n1n1r, nnknnkmaxr);
        double *tmp_n11 = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(tmp_n11, nnkmax);
        double *tmp_n1r = (double *) R_chk_calloc(nnkmaxr, sizeof(double)); zeros(tmp_n1r, nnkmaxr);
        double *tmp_nnknnkmax = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(tmp_nnknnkmax, nnknnkmax);
        double *tmp_nknkmax = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(tmp_nknkmax, nknkmax);
        double *tmp_nkmaxr = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(tmp_nkmaxr, nkmaxr);

        // Set-up storage for sampling for leave-one-out model fit
        double *cv_v_eta = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_v_eta, nnkmax);
        double *cv_v_xi = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(cv_v_xi, nnkmax);
        double *cv_v_beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(cv_v_beta, p);
        double *cv_v_z = (double *) R_chk_calloc(nnkmaxr, sizeof(double)); zeros(cv_v_z, nnkmaxr);

        int cv_index = 0;
        int start_index = 0;
        int end_index = 0;
        int cv_i = 0;
        int sMC_CV = 0;
        int loopd_nMC_nkmax = loopd_nMC * nkmax;
        double *loopd_val_MC_CV = (double *) R_chk_calloc(loopd_nMC_nkmax, sizeof(double)); zeros(loopd_val_MC_CV, loopd_nMC_nkmax);

        GetRNGstate();

        for(cv_index = 0; cv_index < CV_K; cv_index++){

          // set-up partition sizes and indices
          nk = sizesCV[cv_index];
          nknk = nk * nk;
          nnk = n - nk;
          nnknnk = nnk * nnk;
          nnkr = nnk * r;
          nkr = nk * r;

          start_index = startsCV[cv_index];
          end_index = endsCV[cv_index];
          // Rprintf("CV index: %d, start index: %d, end index: %d\n", cv_index, start_index, end_index);

          // Block-deleted data
          copyVecExcludingBlock(Y, cvY, n, start_index, end_index);                                                 // Block-deleted Y
          copyVecExcludingBlock(nBinom, cv_nBinom, n, start_index, end_index);                                      // Block-deleted nBinom
          copyMatrixDelRowBlock(X, n, p, cvX, start_index, end_index);                                              // Block-deleted X
          copyMatrixDelRowBlock(X_tilde, n, r, cvX_tilde, start_index, end_index);                                  // Block-deleted X_tilde

          // Held-out data
          copyMatrixRowBlock(X, n, p, X_pred, start_index, end_index);                                              // Held-out X = X_pred
          copyMatrixRowBlock(X_tilde, n, r, X_tilde_pred, start_index, end_index);                                  // Held-out X_tilde = X_tilde_pred
          copyVecBlock(Y, Y_pred, n, start_index, end_index);                                                       // Held-out Y = Y_pred
          copyVecBlock(nBinom, nBinom_pred, n, start_index, end_index);                                             // Held-out nBinom = nBinom_pred

          // Cross-validated XtX, substract X[i]*t(X[i]) from XtX, instead of multiplying again
          F77_NAME(dgemm)(ytran, ntran, &p, &p, &nk, &one, X_pred, &nk, X_pred, &nk, &zero, cvXtX, &p FCONE FCONE);
          F77_NAME(dscal)(&pp, &negOne, cvXtX, &incOne);
          F77_NAME(daxpy)(&pp, &one, XtX, &incOne, cvXtX, &incOne);

          // Cross-validated XTildetX, delete i-th row block from every n-th block of original XTildetX, instead of multiplying again
          copyMatrixDelRowBlock_vc(XTildetX, nr, p, cvXTildetX, start_index, end_index, n);

          // Constructing cross-validated Vz for each spatial-temporal process model, and also Schur complement for prediction
          if(processType == "independent.shared" || processType == "multivariate"){

            // spatial-temporal covariance matrix
            copyMatrixDelRowColBlock(Vz, n, n, cvVz, start_index, end_index, start_index, end_index);
            cholBlockDelUpdate(n, cholVz, start_index, end_index, cvCholVz, tmp_nnknnkmax, tmp_n11);

            // Pre-processing for spatial prediction
            copyMatrixColDelRowBlock(Vz, n, n, cvCz, start_index, end_index, start_index, end_index);                          // cvCz = Vz[-ids, ids]
            F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, cvCholVz, &nnk, cvCz, &nnk FCONE FCONE FCONE FCONE);  // cvCz <- inv(Lz)*cvCz
            F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, cvCz, &nnk, cvCz, &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE);  // tmp_nknkmax = t(Cz)*inv(Vz)*Cz
            copyMatrixRowColBlock(Vz, n, n, z_tilde_cov, start_index, end_index, start_index, end_index);                      // z_tilde_cov = Vz[ids, ids]
            F77_NAME(daxpy)(&nknk, &negOne, tmp_nknkmax, &incOne, z_tilde_cov, &incOne);
            F77_NAME(dpotrf)(lower, &nk, z_tilde_cov, &nk, &info FCONE); if(info != 0){perror("c++ error: z_schur dpotrf failed\n");}
            mkLT(z_tilde_cov, nk);

          }else if(processType == "independent"){

            for(k = 0; k < r; k++){

              // spatial-temporal covariance matrix
              copyMatrixDelRowColBlock(&Vz[nn * k], n, n, &cvVz[nnknnk * k], start_index, end_index, start_index, end_index);
              cholBlockDelUpdate(n, &cholVz[nn * k], start_index, end_index, &cvCholVz[nnknnk * k], tmp_nnknnkmax, tmp_n11);

              // Pre-processing for spatial prediction
              copyMatrixColDelRowBlock(&Vz[nn * k], n, n, &cvCz[nnk * nk * k], start_index, end_index, start_index, end_index);
              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, &cvCholVz[nnknnk * k], &nnk, &cvCz[nnk * nk * k], &nnk FCONE FCONE FCONE FCONE);
              F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, &cvCz[nnk * nk * k], &nnk, &cvCz[nnk * nk * k], &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE);
              copyMatrixRowColBlock(&Vz[nn * k], n, n, &z_tilde_cov[nknk * k], start_index, end_index, start_index, end_index);
              F77_NAME(daxpy)(&nknk, &negOne, tmp_nknkmax, &incOne, &z_tilde_cov[nknk * k], &incOne);
              F77_NAME(dpotrf)(lower, &nk, &z_tilde_cov[nknk * k], &nk, &info FCONE); if(info != 0){perror("c++ error: z_schur dpotrf failed\n");}
              mkLT(&z_tilde_cov[nknk * k], nk);

            }

          }else if(processType == "multivariate2"){

            // spatial-temporal covariance matrix
            copyMatrixDelRowColBlock_vc(Vz, nr, nr, cvVz, start_index, end_index, start_index, end_index, n);
            cholBlockDelUpdate(n, cholR, start_index, end_index, cvCholR, tmp_nnknnkmax, tmp_n11);
            chol_kron(r, nnk, chol_iwScale, cvCholR, cvCholVz);
            mkLT(cvCholVz, nnkr);

            // Pre-processing for spatial prediction
            copyMatrixColDelRowBlock(R, n, n, cvCz, start_index, end_index, start_index, end_index);
            F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, cvCholR, &nnk, cvCz, &nnk FCONE FCONE FCONE FCONE);
            F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, cvCz, &nnk, cvCz, &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE);
            copyMatrixRowColBlock(R, n, n, z_tilde_cov, start_index, end_index, start_index, end_index);
            F77_NAME(daxpy)(&nknk, &negOne, tmp_nknkmax, &incOne, z_tilde_cov, &incOne);
            F77_NAME(dpotrf)(lower, &nk, z_tilde_cov, &nk, &info FCONE); if(info != 0){perror("c++ error: z_schur dpotrf failed\n");}
            mkLT(z_tilde_cov, nk);

          }

          // Constructing leave-one-out Cholesky factor for I + XTilde*Vz*t(XTilde)
          // It can be shown that it is equivalent to updating Cholesky factor
          // after removing i-th row and i-th column blocks from original chol(I + XTilde*Vz*t(XTilde))
          cholBlockDelUpdate(n, cholIplusXTildeVzXTildet, start_index, end_index, cvCholIplusXTildeVzXTildet, tmp_nnknnkmax, tmp_n11);

          primingGLMvc(nnk, p, r, cvX, cvX_tilde, cvXtX, cvXTildetX, VbetaInv, cvVz, processType, cvCholIplusXTildeVzXTildet,
                       sigmaSq_xi, tmp_n1n1r, cvD1Inv, cvD1InvB1, cvCholschurA1, cvDInvB_pn, cvDInvB_nrn, cvCholschurA);

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
              cv_v_beta[j] = rnorm(0.0, dtemp3);
            }

            for(cv_i = 0; cv_i < nnk; cv_i++){
              cv_v_xi[cv_i] = rnorm(0.0, sigma_xi);
            }

            if(processType == "independent.shared"){
              dtemp1 = 0.5 * nu_z;
              dtemp2 = 1.0 / dtemp1;
              dtemp3 = rgamma(dtemp1, dtemp2);
              dtemp3 = 1.0 / dtemp3;
              dtemp3 = sqrt(dtemp3);
              for(k = 0; k < r; k++){
                for(cv_i = 0; cv_i < nnk; cv_i++){
                  cv_v_z[k*nnk + cv_i] = rnorm(0.0, dtemp3);
                }
              }
            }else if(processType == "independent"){
              for(k = 0; k < r; k++){
                dtemp1 = 0.5 * nu_z;
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp3 = 1.0 / dtemp3;
                dtemp3 = sqrt(dtemp3);
                for(cv_i = 0; cv_i < nnk; cv_i++){
                  cv_v_z[k*nnk + cv_i] = rnorm(0.0, dtemp3);
                }
              }
            }else if(processType == "multivariate"){

              for(k = 0; k < r; k++){
                for(cv_i = 0; cv_i < nnk; cv_i++){
                  tmp_n1r[k*nnk + cv_i] = rnorm(0.0, 1.0);
                }
              }
              rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr);
              F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){perror("c++ error: samp_Sigma dpotrf failed\n");}
              mkLT(samp_Sigma, r);
              F77_NAME(dgemm)(ntran, ytran, &nnk, &r, &r, &one, tmp_n1r, &nnk, samp_Sigma, &r, &zero, cv_v_z, &nnk FCONE FCONE);

            }else if(processType == "multivariate2"){
              for(k = 0; k < r; k++){
                dtemp1 = 0.5 * nu_z;
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp3 = 1.0 / dtemp3;
                dtemp3 = sqrt(dtemp3);
                for(cv_i = 0; cv_i < nnk; cv_i++){
                  cv_v_z[k*nnk + cv_i] = rnorm(0.0, dtemp3);
                }
              }
            }

            // projection step
            projGLMvc(nnk, p, r, cvX, cvX_tilde, sigmaSq_xi, Lbeta, cvCholVz, processType,
                      cv_v_eta, cv_v_xi, cv_v_beta, cv_v_z, cvD1Inv, cvD1InvB1, cvCholschurA1,
                      cvDInvB_pn, cvDInvB_nrn, cvCholschurA, tmp_n1r);

            // Prediction at held-out point for each spatial-temporal process model
            if(processType == "independent.shared"){

              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &r, &one, cvCholVz, &nnk, cv_v_z, &nnk FCONE FCONE FCONE FCONE);  // cv_v_z <- inv(Lz)*v_z
              F77_NAME(dgemm)(ytran, ntran, &nk, &r, &nnk, &one, cvCz, &nnk, cv_v_z, &nnk, &zero, z_tilde_mu, &nk FCONE FCONE);   // z_tilde_mu <- t(Cz)*inv(Vz)*v_z
              for(k = 0; k < r; k++){
                PCM_dist[0] = pow(F77_NAME(dnrm2)(&nnk, &cv_v_z[nnk * k], &incOne), 2);

                // sample z_tilde
                dtemp1 = 0.5 * (nu_z + nnk);
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp2 = 1.0 / dtemp3;
                dtemp1 = (PCM_dist[0] + nu_z) / (nu_z + nnk);
                dtemp3 = dtemp1 * dtemp2;
                dtemp1 = sqrt(dtemp3);
                for(cv_i = 0; cv_i < nk; cv_i++){
                  z_tilde[k*nk + cv_i] = rnorm(0.0, dtemp1);
                }
                F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
                F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, tmp_nkmaxr, &incOne);
                F77_NAME(dcopy)(&nkr, tmp_nkmaxr, &incOne, z_tilde, &incOne);
              }

            }else if(processType == "independent"){

              for(k = 0; k < r; k++){
                F77_NAME(dtrsv)(lower, ntran, nunit, &nnk, &cvCholVz[nnknnk * k], &nnk, &cv_v_z[nnk * k], &incOne FCONE FCONE FCONE);
                F77_NAME(dgemv)(ytran, &nnk, &nk, &one, &cvCz[nnk * nk * k], &nnk, &cv_v_z[nnk * k], &incOne, &zero, &z_tilde_mu[nk * k], &incOne FCONE);
                PCM_dist[0] = pow(F77_NAME(dnrm2)(&nnk, &cv_v_z[nnk * k], &incOne), 2);

                // sample z_tilde
                dtemp1 = 0.5 * (nu_z + nnk);
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp2 = 1.0 / dtemp3;
                dtemp1 = (PCM_dist[0] + nu_z) / (nu_z + nnk);
                dtemp3 = dtemp1 * dtemp2;
                dtemp1 = sqrt(dtemp3);
                for(cv_i = 0; cv_i < nk; cv_i++){
                  z_tilde[k*nk + cv_i] = rnorm(0.0, dtemp1);
                }
                F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
                F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, tmp_nkmaxr, &incOne);
                F77_NAME(dcopy)(&nkr, tmp_nkmaxr, &incOne, z_tilde, &incOne);
              }

            }else if(processType == "multivariate"){

              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &r, &one, cvCholVz, &nnk, cv_v_z, &nnk FCONE FCONE FCONE FCONE);  // cv_v_z <- invchol(R)*v_z
              F77_NAME(dgemm)(ytran, ntran, &nk, &r, &nnk, &one, cvCz, &nnk, cv_v_z, &nnk, &zero, z_tilde_mu, &nk FCONE FCONE);   // z_tilde_mu <- t(C)*inv(R)*v_z
              F77_NAME(dgemm)(ytran, ntran, &r, &r, &nnk, &one, cv_v_z, &nnk, cv_v_z, &nnk, &zero, PCM_dist, &r FCONE FCONE);     // PCM_dist <- t(v_z)*inv(R)*v_z
              F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, PCM_dist, &incOne);                                                    // PCM_dist <- iwScale + t(v_z)*inv(R)*v_z
              F77_NAME(dpotrf)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              F77_NAME(dpotri)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotri failed\n");}
              F77_NAME(dpotrf)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              mkLT(PCM_dist, r);
              rInvWishart(r, nu_z + nnk + 2*r, PCM_dist, samp_Sigma, tmp_rr);
              F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){perror("c++ error: samp_Sigma dpotrf failed\n");}
              mkLT(samp_Sigma, r);

              for(k = 0; k < r; k++){
                for(cv_i = 0; cv_i < nk; cv_i++){
                  z_tilde[k*nk + cv_i] = rnorm(0.0, 1.0);
                }
              }
              F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
              F77_NAME(dgemm)(ntran, ytran, &nk, &r, &r, &one, tmp_nkmaxr, &nk, samp_Sigma, &r, &zero, z_tilde, &nk FCONE FCONE);
              F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

            }else if(processType == "multivariate2"){

              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &r, &one, cvCholR, &nnk, cv_v_z, &nnk FCONE FCONE FCONE FCONE);   // cv_v_z <- invchol(R)*v_z
              F77_NAME(dgemm)(ytran, ntran, &nk, &r, &nnk, &one, cvCz, &nnk, cv_v_z, &nnk, &zero, z_tilde_mu, &nk FCONE FCONE);   // z_tilde_mu <- t(C)*inv(R)*v_z
              F77_NAME(dgemm)(ytran, ntran, &r, &r, &nnk, &one, cv_v_z, &nnk, cv_v_z, &nnk, &zero, PCM_dist, &r FCONE FCONE);     // PCM_dist <- t(v_z)*inv(R)*v_z
              F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, PCM_dist, &incOne);                                                    // PCM_dist <- iwScale + t(v_z)*inv(R)*v_z
              F77_NAME(dpotrf)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){perror("c++ error: post_iwScale dpotrf failed\n");}
              mkLT(PCM_dist, r);
              for(k = 0; k < r; k++){
                dtemp1 = 0.5 * (nu_z + nnk);
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp2 = 1.0 / dtemp3;
                dtemp1 = sqrt(dtemp3);
                for(cv_i = 0; cv_i < nk; cv_i++){
                  z_tilde[k*nk + cv_i] = rnorm(0.0, dtemp1);
                }
              }
              F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
              F77_NAME(dgemm)(ntran, ytran, &nk, &r, &r, &one, tmp_nkmaxr, &nk, PCM_dist, &r, &zero, z_tilde, &nk FCONE FCONE);
              F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

            }

            // Find canonical parameter (X*beta + X_tilde*z_tilde)
            lmulm_XTilde_VC(ntran, nk, r, 1, X_tilde_pred, z_tilde, tmp_nkmaxr);
            F77_NAME(dgemv)(ntran, &nk, &p, &one, X_pred, &nk, cv_v_beta, &incOne, &one, tmp_nkmaxr, &incOne FCONE);

            // Find CV-LOO-PD
            if(family == family_poisson){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = exp(tmp_nkmaxr[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dpois(Y_pred[cv_i], dtemp1, 1);
              }
            }

            if(family == family_binomial){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = inverse_logit(tmp_nkmaxr[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dbinom(Y_pred[cv_i], nBinom_pred[cv_i], dtemp1, 1);
              }
            }

            if(family == family_binary){
              for(cv_i = 0; cv_i < nk; cv_i++){
                dtemp1 = inverse_logit(tmp_nkmaxr[cv_i]);
                loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV] = dbinom(Y_pred[cv_i], 1.0, dtemp1, 1);
              }
            }


          }

          for(cv_i = 0; cv_i < nk; cv_i++){
            REAL(loopd_out_r)[start_index + cv_i] = logMeanExp(&loopd_val_MC_CV[cv_i*loopd_nMC], loopd_nMC);
          }

        }


        R_chk_free(startsCV);
        R_chk_free(endsCV);
        R_chk_free(sizesCV);
        R_chk_free(cvY);
        R_chk_free(cv_nBinom);
        R_chk_free(cvX);
        R_chk_free(cvX_tilde);
        R_chk_free(X_pred);
        R_chk_free(X_tilde_pred);
        R_chk_free(Y_pred);
        R_chk_free(nBinom_pred);
        R_chk_free(cvVz);
        R_chk_free(cvCholVz);
        R_chk_free(cvCholR);
        R_chk_free(cvCz);
        R_chk_free(z_tilde_cov);
        R_chk_free(z_tilde_mu);
        R_chk_free(z_tilde);
        R_chk_free(PCM_dist);
        R_chk_free(cvXtX);
        R_chk_free(cvXTildetX);
        R_chk_free(cvCholIplusXTildeVzXTildet);
        R_chk_free(cvD1Inv);
        R_chk_free(cvD1InvB1);
        R_chk_free(cvCholschurA1);
        R_chk_free(cvDInvB_pn);
        R_chk_free(cvDInvB_nrn);
        R_chk_free(cvCholschurA);
        R_chk_free(tmp_n1n1r);
        R_chk_free(tmp_n11);
        R_chk_free(tmp_n1r);
        R_chk_free(tmp_nnknnkmax);
        R_chk_free(tmp_nknkmax);
        R_chk_free(tmp_nkmaxr);
        R_chk_free(cv_v_eta);
        R_chk_free(cv_v_xi);
        R_chk_free(cv_v_beta);
        R_chk_free(cv_v_z);
        R_chk_free(loopd_val_MC_CV);

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

      // samples of z
      SET_VECTOR_ELT(result_r, 2, samples_xi_r);
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("xi"));

      Rf_namesgets(result_r, resultName_r);

    }

    UNPROTECT(nProtect);
    // return R_NilValue;

    return result_r;

  } // end stvcGLMexactLOO

}