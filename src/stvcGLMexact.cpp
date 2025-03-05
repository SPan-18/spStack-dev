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

  SEXP stvcGLMexact(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                    SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                    SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                    SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                    SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, k, s, info, nProtect = 0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    const double one = 1.0;
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

        if(processType == "independent.shared"){

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

        }else if(processType == "multivariate"){

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
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                     // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);                                                     // Lbeta = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");} // VbetaInv = chol2inv(Vbeta)

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

      for(k = 0; k < r; k++){
        dtemp1 = 0.5 * nu_z;
        dtemp2 = 1.0 / dtemp1;
        dtemp3 = rgamma(dtemp1, dtemp2);
        dtemp3 = 1.0 / dtemp3;
        dtemp3 = sqrt(dtemp3);
        for(i = 0; i < n; i++){
            v_z[k*n + i] = rnorm(0.0, dtemp3);                                           // v_z ~ t
        }
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

    R_chk_free(tmp_nr);

    R_chk_free(D1Inv);
    R_chk_free(D1InvB1);
    R_chk_free(cholschurA1);
    R_chk_free(DInvB_pn);
    R_chk_free(DInvB_nrn);
    R_chk_free(cholschurA);

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

  } // end stvcGLMexact
}