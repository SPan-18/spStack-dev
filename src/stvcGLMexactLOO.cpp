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

// Fits of the candidate models sharing the process parameters (phi_s, phi_t) that differ only in the boundary
// adjustment parameter epsilon (a vector epsilon_r of nEps values), with optional leave-one-out predictive densities.
// The pre-processing of the full data and of each leave-one-out or cross-validation subset does not depend on
// epsilon, so it is computed once and shared by the nEps fits: the posterior samples of the nEps models are drawn
// first (in the order of epsilon_r), and then, for each held-out site or fold, the Monte Carlo draws of the nEps
// models are made in turn. With nEps = 1 the random-number stream is that of a single fit.
// cvUpdate (K-fold CV only): 1 = the pre-processing of each block-deleted data set is obtained from the full-data
// one by deletion updates (scalar loops); 0 = it is recomputed directly (level-3 BLAS; faster with an optimized
// multi-threaded BLAS). Both give the same result up to floating-point rounding.
// Returns a list of nEps fits.
static SEXP stvcGLMexactLOO_fit(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                                SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                                SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                                SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                                SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                                SEXP CV_K_r, SEXP loopd_nMC_r, SEXP cvUpdate_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, k, s, e, info, nProtect = 0;
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

    // supported spatial-temporal process models
    if(processType != "independent.shared" && processType != "independent" && processType != "multivariate"){
      Rf_error("c++ error: process.type must be one of 'independent', 'independent.shared' or 'multivariate'.");
    }
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    double *Vz = NULL;

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

        }
    }

    // boundary adjustment parameters
    int nEps = Rf_length(epsilon_r);
    double *epsVec = REAL(epsilon_r);
    double epsilon = epsVec[0];

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];

    // Leave-one-out predictive density details
    int loopd = INTEGER(loopd_r)[0];
    std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));
    int CV_K = INTEGER(CV_K_r)[0];
    int loopd_nMC = INTEGER(loopd_nMC_r)[0];
    int cvUpdate = INTEGER(cvUpdate_r)[0];

    const char *exact_str = "exact";
    const char *cv_str = "cv";
  
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
      if(nEps == 1){
        Rprintf("\tBoundary adjustment parameter = %.2f.\n\n", epsilon);
      }else{
        Rprintf("\tBoundary adjustment parameters =");
        for(e = 0; e < nEps; e++){
          Rprintf(" %.2f", epsVec[e]);
        }
        Rprintf(".\n\n");
      }

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
      }
      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/

    double *cholVz = NULL;               // define NULL pointer for chol(Vz)
    double *chol_iwScale = NULL;

    // Find Cholesky of Vz
    if(processType == "independent.shared"){

        cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);            // nxn matrix chol(Vz)
        F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i); it is numerically singular, check for nearly coincident space-time locations or very small phi_s, phi_t.\n", info);}
        mkLT(cholVz, n);

    }else if(processType == "independent"){

        cholVz = (double *) R_alloc(nnr, sizeof(double)); zeros(cholVz, nnr);          // r nxn matrices chol(Vz)
        F77_NAME(dcopy)(&nnr, Vz, &incOne, cholVz, &incOne);
        for(k = 0; k < r; k++){
            F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE); if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i); it is numerically singular, check for nearly coincident space-time locations or very small phi_s, phi_t.\n", info);}
            mkLT(&cholVz[nn * k], n);
        }

    }else if(processType == "multivariate"){

      cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);            // nxn matrix chol(Vz)
      F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
      F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){Rf_error("c++ error: Cholesky factorization of the spatial-temporal correlation matrix failed (info = %i); it is numerically singular, check for nearly coincident space-time locations or very small phi_s, phi_t.\n", info);}
      mkLT(cholVz, n);
      chol_iwScale = (double *) R_alloc(rr, sizeof(double)); zeros(chol_iwScale, rr);
      F77_NAME(dcopy)(&rr, iwScale, &incOne, chol_iwScale, &incOne);
      F77_NAME(dpotrf)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){Rf_error("c++ error: the inverse-Wishart scale matrix iw.scale is not positive definite.\n");}
      F77_NAME(dpotri)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){Rf_error("c++ error: the inverse-Wishart scale matrix iw.scale is not positive definite.\n");} // chol_iwScale = chol2inv(iwScale)
      F77_NAME(dpotrf)(lower, &r, chol_iwScale, &r, &info FCONE); if(info != 0){Rf_error("c++ error: the inverse-Wishart scale matrix iw.scale is not positive definite.\n");}
      mkLT(chol_iwScale, r);

    }

    // diagnostics: for each distinct correlation matrix (r of them for 'independent'), the correlations of the
    // farthest-apart and the closest space-time locations and the smallest relative Cholesky pivot of its factor
    // (no extra factorization)
    int nDiag = (processType == "independent") ? r : 1;
    double *diagPivot = (double *) R_alloc(nDiag, sizeof(double));
    double *diagMinCor = (double *) R_alloc(nDiag, sizeof(double));
    double *diagMaxCor = (double *) R_alloc(nDiag, sizeof(double));
    for(k = 0; k < nDiag; k++){
      corOffDiagRange(&Vz[nn * k], n, &diagMinCor[k], &diagMaxCor[k]);
      diagPivot[k] = minRelPivot(&cholVz[nn * k], n, NULL, 1.0);
    }

    // Allocations for VbetaInv
    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);           // allocate VbetaInv
    double *Lbeta = (double *) R_alloc(pp, sizeof(double)); zeros(Lbeta, pp);                 // Cholesky of Vbeta

    // Find VbetaInv
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                           // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){Rf_error("c++ error: prior covariance of beta is not positive definite.\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dcopy)(&pp, VbetaInv, &incOne, Lbeta, &incOne);                                                           // Lbeta = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){Rf_error("c++ error: inversion of the prior covariance of beta failed.\n");}       // VbetaInv = chol2inv(Vbeta)

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
    if(info != 0){Rf_error("c++ error: Cholesky factorization of I + XTilde*Vz*t(XTilde) failed (info = %i).\n", info);}
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
    info = primingGLMvc(n, p, r, X, X_tilde, VbetaInv, Vz, processType, cholIplusXTildeVzXTildet,
                        sigmaSq_xi, tmp_nnr, D1Inv, D1InvB1, cholschurA1, DInvB_pn, DInvB_nrn, cholschurA);

    R_chk_free(tmp_nnr);
    if(info != 0){
      R_chk_free(D1Inv); R_chk_free(D1InvB1); R_chk_free(cholschurA1); R_chk_free(DInvB_pn);
      R_chk_free(DInvB_nrn); R_chk_free(cholschurA);
      glmPrimingError(info);
    }

    // failure in the sampling, leave-one-out or cross-validation loops: the loop is left (goto), its heap memory is
    // freed, and the error is raised before the return object is made (codes as in glmLOOError)
    int failCode = 0;

    /*****************************************
     Set-up posterior sampling
     *****************************************/
    // posterior samples of sigma-sq and beta
    // (one entry per epsilon)
    SEXP samples_beta_l = PROTECT(Rf_allocVector(VECSXP, nEps)); nProtect++;
    SEXP samples_z_l = PROTECT(Rf_allocVector(VECSXP, nEps)); nProtect++;
    SEXP samples_xi_l = PROTECT(Rf_allocVector(VECSXP, nEps)); nProtect++;
    for(e = 0; e < nEps; e++){
      SET_VECTOR_ELT(samples_beta_l, e, Rf_allocMatrix(REALSXP, p, nSamples));
      SET_VECTOR_ELT(samples_z_l, e, Rf_allocMatrix(REALSXP, nr, nSamples));
      SET_VECTOR_ELT(samples_xi_l, e, Rf_allocMatrix(REALSXP, n, nSamples));
    }

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

    // Posterior samples are drawn in blocks of nBlockMax: within a block the random variates are drawn in exactly
    // the order of a draw-by-draw loop, directly into the output matrices, and the block is then projected at once
    // with level-3 BLAS (projGLMvcbatch).
    const int nBlockMax = 64;
    int nBlock = 0, bb = 0;

    int nnBlockMax = n * nBlockMax;
    int pnBlockMax = p * nBlockMax;
    int nrnBlockMax = nr * nBlockMax;
    double *V_eta = (double *) R_chk_calloc(nnBlockMax, sizeof(double)); zeros(V_eta, nnBlockMax);
    double *tmp_nb = (double *) R_chk_calloc(nnBlockMax, sizeof(double)); zeros(tmp_nb, nnBlockMax);
    double *tmp_pb = (double *) R_chk_calloc(pnBlockMax, sizeof(double)); zeros(tmp_pb, pnBlockMax);
    double *tmp_nrb = (double *) R_chk_calloc(nrnBlockMax, sizeof(double)); zeros(tmp_nrb, nrnBlockMax);
    double *V_beta = NULL, *V_xi = NULL, *V_z = NULL;

    // posterior samples of the nEps models, one after the other
    for(e = 0; e < nEps; e++){

      epsilon = epsVec[e];

      for(s = 0; s < nSamples; s += nBlockMax){

        nBlock = std::min(nBlockMax, nSamples - s);
        V_beta = &REAL(VECTOR_ELT(samples_beta_l, e))[(R_xlen_t) s * p];
        V_xi = &REAL(VECTOR_ELT(samples_xi_l, e))[(R_xlen_t) s * n];
        V_z = &REAL(VECTOR_ELT(samples_z_l, e))[(R_xlen_t) s * nr];

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

            for(i = 0; i < n; i++){
              V_xi[bb*n + i] = rnorm(0.0, sigma_xi);                                                  // v_xi ~ N(0, sigmaSq_xi)
            }

            dtemp1 = 0.5 * nu_beta;
            dtemp2 = 1.0 / dtemp1;
            dtemp3 = rgamma(dtemp1, dtemp2);
            dtemp3 = 1.0 / dtemp3;
            dtemp3 = sqrt(dtemp3);
            for(j = 0; j < p; j++){
              V_beta[bb*p + j] = rnorm(0.0, dtemp3);                                                  // v_beta ~ t
            }

            if(processType == "independent.shared"){
              dtemp1 = 0.5 * nu_z;
              dtemp2 = 1.0 / dtemp1;
              dtemp3 = rgamma(dtemp1, dtemp2);
              dtemp3 = 1.0 / dtemp3;
              dtemp3 = sqrt(dtemp3);
              for(k = 0; k < r; k++){
                for(i = 0; i < n; i++){
                  V_z[bb*nr + k*n + i] = rnorm(0.0, dtemp3);
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
                  V_z[bb*nr + k*n + i] = rnorm(0.0, dtemp3);
                }
              }
            }else if(processType == "multivariate"){

              for(k = 0; k < r; k++){
                for(i = 0; i < n; i++){
                  tmp_nr[k*n + i] = rnorm(0.0, 1.0);
                }
              }
              if(rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr) != 0){ failCode = 5; goto fit_done; }
              F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){ failCode = 5; goto fit_done; }
              mkLT(samp_Sigma, r);                                                                       // zero the upper triangle
              F77_NAME(dgemm)(ntran, ytran, &n, &r, &r, &one, tmp_nr, &n, samp_Sigma, &r, &zero, &V_z[bb*nr], &n FCONE FCONE);

            }

        }

        // projection step for the block
        projGLMvcbatch(n, p, r, nBlock, X, X_tilde, sigmaSq_xi, Lbeta, cholVz, processType,
                       V_eta, V_xi, V_beta, V_z, D1Inv, D1InvB1, cholschurA1,
                       DInvB_pn, DInvB_nrn, cholschurA, tmp_nrb, tmp_nb, tmp_pb);

      }

    }

    fit_done:

    R_chk_free(V_eta);
    R_chk_free(tmp_nb);
    R_chk_free(tmp_pb);
    R_chk_free(tmp_nrb);

    PutRNGstate();

    R_chk_free(v_eta);
    R_chk_free(v_xi);
    R_chk_free(v_beta);
    R_chk_free(v_z);
    R_chk_free(tmp_nr);

    R_chk_free(cholschurA1);
    // D1Inv, D1InvB1, DInvB_pn, DInvB_nrn and cholschurA are kept: the exact LOO and CV pre-processing is
    // obtained from them by deletion updates (cholSchurGLMvcDel); they are freed at the end

    // make return object
    SEXP result_r, resultName_r;
    SEXP loopd_out_l = R_NilValue;                                                       // leave-one-out predictive densities, one per epsilon

    if(loopd && failCode == 0){

      if(verbose){
        Rprintf("Evaluating leave-one-out predictive densities.\n");
      }

      loopd_out_l = PROTECT(Rf_allocVector(VECSXP, nEps)); nProtect++;
      for(e = 0; e < nEps; e++){
        SET_VECTOR_ELT(loopd_out_l, e, Rf_allocVector(REALSXP, n));
      }

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
        double *looCholVz = NULL;
        double *looCz = NULL;

        if(corfn == "gneiting-decay"){

          if(processType == "independent.shared" || processType == "multivariate"){

            looCholVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVz, n1n1);
            looCz = (double *) R_chk_calloc(n1, sizeof(double)); zeros(looCz, n1);

          }else if(processType == "independent"){

            looCholVz = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(looCholVz, n1n1r);
            looCz = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(looCz, n1r);

          }

        }


        // set-up pre-processing memory allocations for priming on leave-one-out data
        double *looD1Inv = (double *) R_chk_calloc(n1rn1r, sizeof(double)); zeros(looD1Inv, n1rn1r);
        double *looD1InvB1 = (double *) R_chk_calloc(n1rp, sizeof(double)); zeros(looD1InvB1, n1rp);
        double *looCholschurA1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(looCholschurA1, pp);
        double *looDInvB_pn = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looDInvB_pn, n1p);
        double *looDInvB_nrn = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(looDInvB_nrn, n1n1r);
        double *looCholschurA = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholschurA, n1n1);
        double *tmp_n11 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n11, n1);
        double *tmp_n1r = (double *) R_chk_calloc(n1r, sizeof(double)); zeros(tmp_n1r, n1r);

        // Workspace for the deletion update of the pre-processing (cholSchurGLMvcDel with a block of size 1)
        double *del_PB = (double *) R_chk_calloc(n, sizeof(double)); zeros(del_PB, n);
        double *del_QB = (double *) R_chk_calloc(n, sizeof(double)); zeros(del_QB, n);
        double *del_QBK = (double *) R_chk_calloc(n, sizeof(double)); zeros(del_QBK, n);
        double *del_LP = (double *) R_chk_calloc(1, sizeof(double)); zeros(del_LP, 1);
        double *del_LQ = (double *) R_chk_calloc(1, sizeof(double)); zeros(del_LQ, 1);
        double *del_WB = (double *) R_chk_calloc(p, sizeof(double)); zeros(del_WB, p);
        double *del_Z = (double *) R_chk_calloc(p, sizeof(double)); zeros(del_Z, p);
        double *del_H = (double *) R_chk_calloc(nr, sizeof(double)); zeros(del_H, nr);
        double *del_HK = (double *) R_chk_calloc(nr, sizeof(double)); zeros(del_HK, nr);
        double *del_A2 = (double *) R_chk_calloc(nr, sizeof(double)); zeros(del_A2, nr);
        double *del_tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(del_tmp_np, np);
        double *del_w = (double *) R_chk_calloc(n, sizeof(double)); zeros(del_w, n);

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

        // blocks of nBlockMax Monte Carlo draws
        int n1BlockMax = n1 * nBlockMax;
        int n1rBlockMax = n1r * nBlockMax;
        int pBlockMax = p * nBlockMax;
        int rBlockMax = r * nBlockMax;
        int rrBlockMax = rr * nBlockMax;
        double *LV_eta = (double *) R_chk_calloc(n1BlockMax, sizeof(double)); zeros(LV_eta, n1BlockMax);
        double *LV_xi = (double *) R_chk_calloc(n1BlockMax, sizeof(double)); zeros(LV_xi, n1BlockMax);
        double *LV_tmpn = (double *) R_chk_calloc(n1BlockMax, sizeof(double)); zeros(LV_tmpn, n1BlockMax);
        double *LV_z = (double *) R_chk_calloc(n1rBlockMax, sizeof(double)); zeros(LV_z, n1rBlockMax);
        double *LV_tmpnr = (double *) R_chk_calloc(n1rBlockMax, sizeof(double)); zeros(LV_tmpnr, n1rBlockMax);
        double *LV_beta = (double *) R_chk_calloc(pBlockMax, sizeof(double)); zeros(LV_beta, pBlockMax);
        double *LV_tmpp = (double *) R_chk_calloc(pBlockMax, sizeof(double)); zeros(LV_tmpp, pBlockMax);
        double *LV_gam = (double *) R_chk_calloc(rBlockMax, sizeof(double)); zeros(LV_gam, rBlockMax);
        double *LV_nrm = (double *) R_chk_calloc(rBlockMax, sizeof(double)); zeros(LV_nrm, rBlockMax);
        double *LV_bart = (double *) R_chk_calloc(rrBlockMax, sizeof(double)); zeros(LV_bart, rrBlockMax);

        GetRNGstate();

        for(loo_index = 0; loo_index < n; loo_index++){

          if(pendingInterrupt()){ failCode = 6; goto loo_done; }            // user interrupt: free memory, then stop

          // Prepare leave-one-out data
          copyVecExcludingOne(Y, looY, n, loo_index);                           // Leave-one-out Y
          copyVecExcludingOne(nBinom, loo_nBinom, n, loo_index);                // Leave-one-out nBinom
          copyMatrixDelRow(X, n, p, looX, loo_index);                           // Row-deleted X
          copyMatrixDelRow(X_tilde, n, r, looX_tilde, loo_index);               // Row-deleted X_tilde
          copyMatrixRowToVec(X, n, p, X_pred, loo_index);                       // Copy left out X into X_pred
          copyMatrixRowToVec(X_tilde, n, r, X_tilde_pred, loo_index);           // Copy left out X into X_pred

          // Constructing leave-one-out Vz for each spatial-temporal process model, and also Schur complement for prediction
          if(processType == "independent.shared" || processType == "multivariate"){

            cholRowDelUpdate(n, cholVz, loo_index, looCholVz, tmp_n11);

            copyVecExcludingOne(&Vz[loo_index*n], looCz, n, loo_index);                                            // looCz = Vz[-i,i]
            F77_NAME(dtrsv)(lower, ntran, nunit, &n1, looCholVz, &n1, looCz, &incOne FCONE FCONE FCONE);           // looCz = LzInv * Cz
            dtemp1 = pow(F77_NAME(dnrm2)(&n1, looCz, &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
            z_tilde_var[0] = Vz[loo_index*n + loo_index] - dtemp1;                                                 // z_tilde_var = Vz_tilde - Czt*VzInv*Cz

          }else if(processType == "independent"){

            for(k = 0; k < r; k++){
              cholRowDelUpdate(n, &cholVz[nn*k], loo_index, &looCholVz[n1n1*k], tmp_n11);

              copyVecExcludingOne(&Vz[nn*k + loo_index*n], &looCz[n1*k], n, loo_index);                                     // looCz = Vz[-i,i]
              F77_NAME(dtrsv)(lower, ntran, nunit, &n1, &looCholVz[n1n1*k], &n1, &looCz[n1*k], &incOne FCONE FCONE FCONE);  // looCz = LzInv * Cz
              dtemp1 = pow(F77_NAME(dnrm2)(&n1, &looCz[n1*k], &incOne), 2);                                                 // dtemp1 = Czt*VzInv*Cz
              z_tilde_var[k] = Vz[nn*k + loo_index*n + loo_index] - dtemp1;

            }

          }

          // Pre-processing for projGLMvc() on the leave-one-out data by deletion updates of the full-data outputs
          // (O(n^2 r^2)); chol(I/sigmaSqxi + Q) is row-deleted first and then downdated in place by cholSchurGLMvcDel
          cholRowDelUpdate(n, cholschurA, loo_index, looCholschurA, tmp_n11);
          if(cholSchurGLMvcDel(n, p, r, loo_index, loo_index, X, X_tilde, cholIplusXTildeVzXTildet,
                               D1Inv, D1InvB1, DInvB_pn, DInvB_nrn, VbetaInv,
                               looD1Inv, looD1InvB1, looCholschurA1, looDInvB_pn, looDInvB_nrn, looCholschurA,
                               del_PB, del_QB, del_QBK, del_LP, del_LQ, del_WB, del_Z, del_H, del_HK, del_A2,
                               del_tmp_np, del_w) != 0){

            // not numerically positive definite: recompute directly on the leave-one-out data (O(n^3 r^2)).
            // The leave-one-out Cholesky factor of I + XTilde*Vz*t(XTilde) is the row-deletion update of the
            // full-data factor.
            int looVz_size = n1n1;
            if(processType == "independent"){
              looVz_size = n1n1r;
            }
            double *looVz = (double *) R_chk_calloc(looVz_size, sizeof(double)); zeros(looVz, looVz_size);
            double *looCholIplusXTildeVzXTildet = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholIplusXTildeVzXTildet, n1n1);
            double *tmp_n1n1r = (double *) R_chk_calloc(n1n1r, sizeof(double)); zeros(tmp_n1n1r, n1n1r);

            if(processType == "independent.shared" || processType == "multivariate"){
              copyMatrixDelRowCol(Vz, n, n, looVz, loo_index, loo_index);
            }else if(processType == "independent"){
              for(k = 0; k < r; k++){
                copyMatrixDelRowCol(&Vz[nn*k], n, n, &looVz[n1n1*k], loo_index, loo_index);
              }
            }
            cholRowDelUpdate(n, cholIplusXTildeVzXTildet, loo_index, looCholIplusXTildeVzXTildet, tmp_n11);

            failCode = primingGLMvc(n1, p, r, looX, looX_tilde, VbetaInv, looVz, processType, looCholIplusXTildeVzXTildet,
                                    sigmaSq_xi, tmp_n1n1r, looD1Inv, looD1InvB1, looCholschurA1, looDInvB_pn, looDInvB_nrn, looCholschurA);

            R_chk_free(looVz);
            R_chk_free(looCholIplusXTildeVzXTildet);
            R_chk_free(tmp_n1n1r);
            if(failCode != 0){ goto loo_done; }

          }

          // Monte Carlo LOO-PD of the nEps models, sharing the pre-processing above
          for(e = 0; e < nEps; e++){

            epsilon = epsVec[e];

            int nBlockr = 0;
            double *Zb = NULL;

            for(sMC = 0; sMC < loopd_nMC; sMC += nBlockMax){

              nBlock = std::min(nBlockMax, loopd_nMC - sMC);

              // random variates of the block, in the order of a draw-by-draw loop
              for(bb = 0; bb < nBlock; bb++){


                if(family == family_poisson){
                  for(loo_i = 0; loo_i < n1; loo_i++){
                    dtemp1 = looY[loo_i] + epsilon;
                    dtemp2 = 1.0;
                    LV_eta[bb*n1 + loo_i] = rlogGamma(dtemp1);                              // log(Gamma(y + epsilon, 1)), underflow-safe
                  }
                }

                if(family == family_binomial){
                  for(loo_i = 0; loo_i < n1; loo_i++){
                    dtemp1 = looY[loo_i] + epsilon;
                    dtemp2 = loo_nBinom[loo_i];
                    dtemp2 += 2.0 * epsilon;
                    dtemp2 -= dtemp1;
                    LV_eta[bb*n1 + loo_i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
                  }
                }

                if(family == family_binary){
                  for(loo_i = 0; loo_i < n1; loo_i++){
                    dtemp1 = looY[loo_i] + epsilon;
                    dtemp2 = loo_nBinom[loo_i];
                    dtemp2 += 2.0 * epsilon;
                    dtemp2 -= dtemp1;
                    LV_eta[bb*n1 + loo_i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
                  }
                }

                for(loo_i = 0; loo_i < n1; loo_i++){
                  LV_xi[bb*n1 + loo_i] = rnorm(0.0, sigma_xi);
                }

                dtemp1 = 0.5 * nu_beta;
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp1 = 1.0 / dtemp3;
                dtemp2 = sqrt(dtemp1);
                for(j = 0; j < p; j++){
                  LV_beta[bb*p + j] = rnorm(0.0, dtemp2);                                                  // loo_v_beta ~ t_nu_beta(0, 1)
                }

                if(processType == "independent.shared"){
                  dtemp1 = 0.5 * nu_z;
                  dtemp2 = 1.0 / dtemp1;
                  dtemp3 = rgamma(dtemp1, dtemp2);
                  dtemp3 = 1.0 / dtemp3;
                  dtemp3 = sqrt(dtemp3);
                  for(k = 0; k < r; k++){
                    for(loo_i = 0; loo_i < n1; loo_i++){
                      LV_z[bb*n1r + k*n1 + loo_i] = rnorm(0.0, dtemp3);
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
                      LV_z[bb*n1r + k*n1 + loo_i] = rnorm(0.0, dtemp3);
                    }
                  }
                }else if(processType == "multivariate"){
                  for(k = 0; k < r; k++){
                    for(loo_i = 0; loo_i < n1; loo_i++){
                      tmp_n1r[k*n1 + loo_i] = rnorm(0.0, 1.0);
                    }
                  }
                  if(rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr) != 0){ failCode = 5; goto loo_done; }
                  F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){ failCode = 5; goto loo_done; }
                  mkLT(samp_Sigma, r);
                  F77_NAME(dgemm)(ntran, ytran, &n1, &r, &r, &one, tmp_n1r, &n1, samp_Sigma, &r, &zero, &LV_z[bb*n1r], &n1 FCONE FCONE);

                }

                // variates of the z_tilde draw, scaled below once the projection gives the scales:
                // rgamma(0.5*(nu_z + n1), .) and the standard normal of rnorm(0, sd) = sd*norm_rand() per process
                // (independent processes), or the Bartlett factor of the inverse-Wishart draw and r standard
                // normals (multivariate)
                if(processType == "independent.shared" || processType == "independent"){
                  for(k = 0; k < r; k++){
                    dtemp1 = 0.5 * (nu_z + n1);
                    dtemp2 = 1.0 / dtemp1;
                    LV_gam[bb*r + k] = rgamma(dtemp1, dtemp2);
                    LV_nrm[bb*r + k] = norm_rand();
                  }
                }else if(processType == "multivariate"){
                  rWishartBartlett(r, nu_z + n1 + 2*r, &LV_bart[bb*rr]);
                  for(k = 0; k < r; k++){
                    LV_nrm[bb*r + k] = norm_rand();
                  }
                }

              }

              // projection step for the block
              projGLMvcbatch(n1, p, r, nBlock, looX, looX_tilde, sigmaSq_xi, Lbeta, looCholVz, processType,
                             LV_eta, LV_xi, LV_beta, LV_z, looD1Inv, looD1InvB1, looCholschurA1,
                             looDInvB_pn, looDInvB_nrn, looCholschurA, LV_tmpnr, LV_tmpn, LV_tmpp);

              // LV_z <- inv(Lz)*LV_z for every process block of every draw
              nBlockr = nBlock * r;
              if(processType == "independent.shared" || processType == "multivariate"){
                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n1, &nBlockr, &one, looCholVz, &n1, LV_z, &n1 FCONE FCONE FCONE FCONE);
              }else if(processType == "independent"){
                for(k = 0; k < r; k++){
                  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n1, &nBlock, &one, &looCholVz[n1n1*k], &n1, &LV_z[n1*k], &n1r FCONE FCONE FCONE FCONE);
                }
              }

              // Prediction at held-out point for each draw of the block
              for(bb = 0; bb < nBlock; bb++){

                Zb = &LV_z[bb*n1r];

                if(processType == "independent.shared" || processType == "independent"){

                  for(k = 0; k < r; k++){
                    if(processType == "independent.shared"){
                      z_tilde_mu[k] = F77_CALL(ddot)(&n1, looCz, &incOne, &Zb[n1*k], &incOne);                        // z_tilde_mu = Czt*VzInv*v_z
                    }else{
                      z_tilde_mu[k] = F77_CALL(ddot)(&n1, &looCz[n1*k], &incOne, &Zb[n1*k], &incOne);
                    }
                    Mdist_r[k] = pow(F77_NAME(dnrm2)(&n1, &Zb[n1*k], &incOne), 2);                                    // Mdist = v_zt*VzInv*v_z

                    // sample z_tilde
                    dtemp1 = 1.0 / LV_gam[bb*r + k];
                    dtemp2 = dtemp1 * (Mdist_r[k] + nu_z) / (nu_z + n1);
                    dtemp3 = sqrt(dtemp2);
                    z_tilde[k] = dtemp3 * LV_nrm[bb*r + k];                                                          // = rnorm(0.0, dtemp3)
                    if(processType == "independent.shared"){
                      z_tilde[k] = z_tilde[k] * sqrt(z_tilde_var[0]);
                    }else{
                      z_tilde[k] = z_tilde[k] * sqrt(z_tilde_var[k]);
                    }
                    z_tilde[k] = z_tilde[k] + z_tilde_mu[k];
                  }

                }else if(processType == "multivariate"){

                  F77_NAME(dgemm)(ytran, ntran, &incOne, &r, &n1, &one, looCz, &n1, Zb, &n1, &zero, z_tilde_mu, &incOne FCONE FCONE);   // z_tilde_mu = t(C)*inv(R)*Z
                  F77_NAME(dgemm)(ytran, ntran, &r, &r, &n1, &one, Zb, &n1, Zb, &n1, &zero, Mdist_rr, &r FCONE FCONE);                  // Mdist = t(Z)*inv(R)*Z
                  F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, Mdist_rr, &incOne);                                                      // Mdist = iwScale + t(Z)*inv(R)*Z
                  F77_NAME(dpotrf)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){ failCode = 5; goto loo_done; }
                  F77_NAME(dpotri)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){ failCode = 5; goto loo_done; }
                  F77_NAME(dpotrf)(lower, &r, Mdist_rr, &r, &info FCONE); if(info != 0){ failCode = 5; goto loo_done; }
                  mkLT(Mdist_rr, r);
                  if(invWishartFromBartlett(r, &LV_bart[bb*rr], Mdist_rr, samp_Sigma, tmp_rr) != 0){ failCode = 5; goto loo_done; }  // = rInvWishart(r, nu_z + n1 + 2r, Mdist_rr, ...)
                  F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){ failCode = 5; goto loo_done; }
                  mkLT(samp_Sigma, r);

                  dtemp1 = sqrt(z_tilde_var[0]);
                  for(k = 0; k < r; k++){
                    tmp_r[k] = dtemp1 * LV_nrm[bb*r + k];                                                            // = rnorm(0.0, dtemp1)
                  }
                  F77_NAME(dgemv)(ntran, &r, &r, &one, samp_Sigma, &r, tmp_r, &incOne, &zero, z_tilde, &incOne FCONE);
                  F77_NAME(daxpy)(&r, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

                }

                dtemp1 = F77_CALL(ddot)(&p, X_pred, &incOne, &LV_beta[bb*p], &incOne);
                dtemp1 += F77_CALL(ddot)(&r, X_tilde_pred, &incOne, z_tilde, &incOne);

                // Find predictive densities from canonical parameter dtemp2 = (X*beta + z)
                if(family == family_poisson){
                  dtemp2 = exp(dtemp1);
                  loopd_val_MC[sMC + bb] = dpois(Y[loo_index], dtemp2, 1);
                }

                if(family == family_binomial){
                  dtemp2 = inverse_logit(dtemp1);
                  loopd_val_MC[sMC + bb] = dbinom(Y[loo_index], nBinom[loo_index], dtemp2, 1);
                }

                if(family == family_binary){
                  dtemp2 = inverse_logit(dtemp1);
                  loopd_val_MC[sMC + bb] = dbinom(Y[loo_index], 1.0, dtemp2, 1);
                }

              }

            }

            REAL(VECTOR_ELT(loopd_out_l, e))[loo_index] = logMeanExp(loopd_val_MC, loopd_nMC);

          }

        }

        loo_done:

        PutRNGstate();

        R_chk_free(looY);
        R_chk_free(loo_nBinom);
        R_chk_free(looX);
        R_chk_free(looX_tilde);
        R_chk_free(X_pred);
        R_chk_free(X_tilde_pred);
        R_chk_free(looCholVz);
        R_chk_free(looCz);
        R_chk_free(looD1Inv);
        R_chk_free(looD1InvB1);
        R_chk_free(looCholschurA1);
        R_chk_free(looDInvB_pn);
        R_chk_free(looDInvB_nrn);
        R_chk_free(looCholschurA);
        R_chk_free(del_PB);
        R_chk_free(del_QB);
        R_chk_free(del_QBK);
        R_chk_free(del_LP);
        R_chk_free(del_LQ);
        R_chk_free(del_WB);
        R_chk_free(del_Z);
        R_chk_free(del_H);
        R_chk_free(del_HK);
        R_chk_free(del_A2);
        R_chk_free(del_tmp_np);
        R_chk_free(del_w);
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
        R_chk_free(LV_eta);
        R_chk_free(LV_xi);
        R_chk_free(LV_tmpn);
        R_chk_free(LV_z);
        R_chk_free(LV_tmpnr);
        R_chk_free(LV_beta);
        R_chk_free(LV_tmpp);
        R_chk_free(LV_gam);
        R_chk_free(LV_nrm);
        R_chk_free(LV_bart);

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
        double *cvX_tilde = (double *) R_chk_calloc(nnkmaxr, sizeof(double)); zeros(cvX_tilde, nnkmaxr);     // Store block-deleted X_tilde
        double *X_pred = (double *) R_chk_calloc(nkmaxp, sizeof(double)); zeros(X_pred, nkmaxp);             // Store held-out X
        double *X_tilde_pred = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(X_tilde_pred, nkmaxr); // Store held-out X
        double *Y_pred = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(Y_pred, nkmax);               // Store held-out Y
        double *nBinom_pred = (double *) R_chk_calloc(nkmax, sizeof(double)); zeros(nBinom_pred, nkmax);     // Store held-out X

        // Set-up storage for pre-processing for cross-validated data
        double *cvCholVz = NULL;
        double *cvCz = NULL;
        double *z_tilde_cov = NULL;
        double *z_tilde_mu = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(z_tilde_mu, nkmaxr);
        double *z_tilde = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(z_tilde, nkmaxr);
        double *PCM_dist = (double *) R_chk_calloc(rr, sizeof(double)); zeros(PCM_dist, rr);

        if(corfn == "gneiting-decay"){

          if(processType == "independent.shared" || processType == "multivariate"){

            cvCholVz = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholVz, nnknnkmax);
            cvCz = (double *) R_chk_calloc(nnkmaxnkmax, sizeof(double)); zeros(cvCz, nnkmaxnkmax);
            z_tilde_cov = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(z_tilde_cov, nknkmax);

          }else if(processType == "independent"){

            cvCholVz = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(cvCholVz, nnknnkmaxr);
            cvCz = (double *) R_chk_calloc(nnkmaxnkmaxr, sizeof(double)); zeros(cvCz, nnkmaxnkmaxr);
            z_tilde_cov = (double *) R_chk_calloc(nknkmaxr, sizeof(double)); zeros(z_tilde_cov, nknkmaxr);

          }

        }


        // set-up pre-processing memory allocations for priming on leave-one-out data
        double *cvD1Inv = (double *) R_chk_calloc(nnkmaxrnnkmaxr, sizeof(double)); zeros(cvD1Inv, nnkmaxrnnkmaxr);
        double *cvD1InvB1 = (double *) R_chk_calloc(nnkmaxrp, sizeof(double)); zeros(cvD1InvB1, nnkmaxrp);
        double *cvCholschurA1 = (double *) R_chk_calloc(pp, sizeof(double)); zeros(cvCholschurA1, pp);
        double *cvDInvB_pn = (double *) R_chk_calloc(nnkmaxp, sizeof(double)); zeros(cvDInvB_pn, nnkmaxp);
        double *cvDInvB_nrn = (double *) R_chk_calloc(nnknnkmaxr, sizeof(double)); zeros(cvDInvB_nrn, nnknnkmaxr);
        double *cvCholschurA = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(cvCholschurA, nnknnkmax);
        double *tmp_n11 = (double *) R_chk_calloc(nnkmax, sizeof(double)); zeros(tmp_n11, nnkmax);
        double *tmp_n1r = (double *) R_chk_calloc(nnkmaxr, sizeof(double)); zeros(tmp_n1r, nnkmaxr);
        double *tmp_nnknnkmax = (double *) R_chk_calloc(nnknnkmax, sizeof(double)); zeros(tmp_nnknnkmax, nnknnkmax);
        double *tmp_nknkmax = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(tmp_nknkmax, nknkmax);
        double *tmp_nkmaxr = (double *) R_chk_calloc(nkmaxr, sizeof(double)); zeros(tmp_nkmaxr, nkmaxr);

        // Workspace for the deletion update of the pre-processing (cholSchurGLMvcDel)
        int nnkmax_del = n * nkmax;
        int nrnkmax_del = nr * nkmax;
        double *del_PB = (double *) R_chk_calloc(nnkmax_del, sizeof(double)); zeros(del_PB, nnkmax_del);        // n x max(nk)
        double *del_QB = (double *) R_chk_calloc(nnkmax_del, sizeof(double)); zeros(del_QB, nnkmax_del);
        double *del_QBK = (double *) R_chk_calloc(nnkmax_del, sizeof(double)); zeros(del_QBK, nnkmax_del);
        double *del_LP = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(del_LP, nknkmax);              // max(nk) x max(nk)
        double *del_LQ = (double *) R_chk_calloc(nknkmax, sizeof(double)); zeros(del_LQ, nknkmax);
        double *del_WB = (double *) R_chk_calloc(nkmaxp, sizeof(double)); zeros(del_WB, nkmaxp);                // max(nk) x p
        double *del_Z = (double *) R_chk_calloc(nkmaxp, sizeof(double)); zeros(del_Z, nkmaxp);
        double *del_H = (double *) R_chk_calloc(nrnkmax_del, sizeof(double)); zeros(del_H, nrnkmax_del);        // nr x max(nk)
        double *del_HK = (double *) R_chk_calloc(nrnkmax_del, sizeof(double)); zeros(del_HK, nrnkmax_del);
        double *del_A2 = (double *) R_chk_calloc(nrnkmax_del, sizeof(double)); zeros(del_A2, nrnkmax_del);
        double *del_tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(del_tmp_np, np);
        double *del_w = (double *) R_chk_calloc(n, sizeof(double)); zeros(del_w, n);

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

        // blocks of nBlockMax Monte Carlo draws
        int nnkBlockMax = nnkmax * nBlockMax;
        int nnkrBlockMax = nnkmaxr * nBlockMax;
        int pBlockMax_cv = p * nBlockMax;
        int rBlockMax_cv = r * nBlockMax;
        int rrBlockMax_cv = rr * nBlockMax;
        int nkrBlockMax = nkmaxr * nBlockMax;
        double *CV_eta = (double *) R_chk_calloc(nnkBlockMax, sizeof(double)); zeros(CV_eta, nnkBlockMax);
        double *CV_xi = (double *) R_chk_calloc(nnkBlockMax, sizeof(double)); zeros(CV_xi, nnkBlockMax);
        double *CV_tmpn = (double *) R_chk_calloc(nnkBlockMax, sizeof(double)); zeros(CV_tmpn, nnkBlockMax);
        double *CV_z = (double *) R_chk_calloc(nnkrBlockMax, sizeof(double)); zeros(CV_z, nnkrBlockMax);
        double *CV_tmpnr = (double *) R_chk_calloc(nnkrBlockMax, sizeof(double)); zeros(CV_tmpnr, nnkrBlockMax);
        double *CV_beta = (double *) R_chk_calloc(pBlockMax_cv, sizeof(double)); zeros(CV_beta, pBlockMax_cv);
        double *CV_tmpp = (double *) R_chk_calloc(pBlockMax_cv, sizeof(double)); zeros(CV_tmpp, pBlockMax_cv);
        double *CV_gam = (double *) R_chk_calloc(rBlockMax_cv, sizeof(double)); zeros(CV_gam, rBlockMax_cv);
        double *CV_nrm = (double *) R_chk_calloc(nkrBlockMax, sizeof(double)); zeros(CV_nrm, nkrBlockMax);
        double *CV_bart = (double *) R_chk_calloc(rrBlockMax_cv, sizeof(double)); zeros(CV_bart, rrBlockMax_cv);

        GetRNGstate();

        for(cv_index = 0; cv_index < CV_K; cv_index++){

          if(pendingInterrupt()){ failCode = 6; goto cv_done; }             // user interrupt: free memory, then stop

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

          // Constructing cross-validated Vz for each spatial-temporal process model, and also Schur complement for prediction
          if(processType == "independent.shared" || processType == "multivariate"){

            // spatial-temporal covariance matrix
            if(cvUpdate){
              cholBlockDelUpdate(n, cholVz, start_index, end_index, cvCholVz, tmp_nnknnkmax, tmp_n11);
            }else{
              copyMatrixDelRowColBlock(Vz, n, n, cvCholVz, start_index, end_index, start_index, end_index);
              F77_NAME(dpotrf)(lower, &nnk, cvCholVz, &nnk, &info FCONE); if(info != 0){ failCode = 3; goto cv_done; }
              mkLT(cvCholVz, nnk);
            }

            // Pre-processing for spatial prediction
            copyMatrixColDelRowBlock(Vz, n, n, cvCz, start_index, end_index, start_index, end_index);                          // cvCz = Vz[-ids, ids]
            F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, cvCholVz, &nnk, cvCz, &nnk FCONE FCONE FCONE FCONE);  // cvCz <- inv(Lz)*cvCz
            F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, cvCz, &nnk, cvCz, &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE);  // tmp_nknkmax = t(Cz)*inv(Vz)*Cz
            copyMatrixRowColBlock(Vz, n, n, z_tilde_cov, start_index, end_index, start_index, end_index);                      // z_tilde_cov = Vz[ids, ids]
            F77_NAME(daxpy)(&nknk, &negOne, tmp_nknkmax, &incOne, z_tilde_cov, &incOne);
            F77_NAME(dpotrf)(lower, &nk, z_tilde_cov, &nk, &info FCONE); if(info != 0){ failCode = 4; goto cv_done; }
            mkLT(z_tilde_cov, nk);

          }else if(processType == "independent"){

            for(k = 0; k < r; k++){

              // spatial-temporal covariance matrix
              if(cvUpdate){
                cholBlockDelUpdate(n, &cholVz[nn * k], start_index, end_index, &cvCholVz[nnknnk * k], tmp_nnknnkmax, tmp_n11);
              }else{
                copyMatrixDelRowColBlock(&Vz[nn * k], n, n, &cvCholVz[nnknnk * k], start_index, end_index, start_index, end_index);
                F77_NAME(dpotrf)(lower, &nnk, &cvCholVz[nnknnk * k], &nnk, &info FCONE); if(info != 0){ failCode = 3; goto cv_done; }
                mkLT(&cvCholVz[nnknnk * k], nnk);
              }

              // Pre-processing for spatial prediction
              copyMatrixColDelRowBlock(&Vz[nn * k], n, n, &cvCz[nnk * nk * k], start_index, end_index, start_index, end_index);
              F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nk, &one, &cvCholVz[nnknnk * k], &nnk, &cvCz[nnk * nk * k], &nnk FCONE FCONE FCONE FCONE);
              F77_NAME(dgemm)(ytran, ntran, &nk, &nk, &nnk, &one, &cvCz[nnk * nk * k], &nnk, &cvCz[nnk * nk * k], &nnk, &zero, tmp_nknkmax, &nk FCONE FCONE);
              copyMatrixRowColBlock(&Vz[nn * k], n, n, &z_tilde_cov[nknk * k], start_index, end_index, start_index, end_index);
              F77_NAME(daxpy)(&nknk, &negOne, tmp_nknkmax, &incOne, &z_tilde_cov[nknk * k], &incOne);
              F77_NAME(dpotrf)(lower, &nk, &z_tilde_cov[nknk * k], &nk, &info FCONE); if(info != 0){ failCode = 4; goto cv_done; }
              mkLT(&z_tilde_cov[nknk * k], nk);

            }

          }

          // Pre-processing for projGLMvc() on the block-deleted data
          int cvDirect = !cvUpdate;
          if(cvUpdate){
            // by deletion updates of the full-data outputs (O(n^2 r^2 nk)); chol(I/sigmaSqxi + Q) is block-deleted
            // first and then downdated in place (nk rank-1 downdates)
            cholBlockDelUpdate(n, cholschurA, start_index, end_index, cvCholschurA, tmp_nnknnkmax, tmp_n11);
            cvDirect = (cholSchurGLMvcDel(n, p, r, start_index, end_index, X, X_tilde, cholIplusXTildeVzXTildet,
                                          D1Inv, D1InvB1, DInvB_pn, DInvB_nrn, VbetaInv,
                                          cvD1Inv, cvD1InvB1, cvCholschurA1, cvDInvB_pn, cvDInvB_nrn, cvCholschurA,
                                          del_PB, del_QB, del_QBK, del_LP, del_LQ, del_WB, del_Z, del_H, del_HK, del_A2,
                                          del_tmp_np, del_w) != 0);
          }
          if(cvDirect){

            // directly on the block-deleted data (O(n^3 r^2)): requested (cvUpdate = 0), or the deletion update was
            // not numerically positive definite. In the latter case the block-deleted Cholesky factor of
            // I + XTilde*Vz*t(XTilde) is the block-deletion update of the full-data factor.
            int cvVz_size = nnknnk;
            if(processType == "independent"){
              cvVz_size = nnknnk * r;
            }
            int nnknnkr = nnknnk * r;
            double *cvVz = (double *) R_chk_calloc(cvVz_size, sizeof(double)); zeros(cvVz, cvVz_size);
            double *cvCholIplusXTildeVzXTildet = (double *) R_chk_calloc(nnknnk, sizeof(double)); zeros(cvCholIplusXTildeVzXTildet, nnknnk);
            double *tmp_n1n1r = (double *) R_chk_calloc(nnknnkr, sizeof(double)); zeros(tmp_n1n1r, nnknnkr);

            if(processType == "independent.shared" || processType == "multivariate"){
              copyMatrixDelRowColBlock(Vz, n, n, cvVz, start_index, end_index, start_index, end_index);
            }else if(processType == "independent"){
              for(k = 0; k < r; k++){
                copyMatrixDelRowColBlock(&Vz[nn * k], n, n, &cvVz[nnknnk * k], start_index, end_index, start_index, end_index);
              }
            }
            if(cvUpdate){
              cholBlockDelUpdate(n, cholIplusXTildeVzXTildet, start_index, end_index, cvCholIplusXTildeVzXTildet, tmp_nnknnkmax, tmp_n11);
            }else{
              rmul_Vz_XTildeT(nnk, r, cvX_tilde, cvVz, tmp_n1n1r, processType);                                          // Vz*t(X_tilde)
              lmulm_XTilde_VC(ntran, nnk, r, nnk, cvX_tilde, tmp_n1n1r, cvCholIplusXTildeVzXTildet);                    // X_tilde*Vz*t(X_tilde)
              for(cv_i = 0; cv_i < nnk; cv_i++){
                cvCholIplusXTildeVzXTildet[cv_i*nnk + cv_i] += 1.0;
              }
              F77_NAME(dpotrf)(lower, &nnk, cvCholIplusXTildeVzXTildet, &nnk, &info FCONE);
              mkLT(cvCholIplusXTildeVzXTildet, nnk);
              if(info != 0){ failCode = 3; }
            }

            if(failCode == 0){
              failCode = primingGLMvc(nnk, p, r, cvX, cvX_tilde, VbetaInv, cvVz, processType, cvCholIplusXTildeVzXTildet,
                                      sigmaSq_xi, tmp_n1n1r, cvD1Inv, cvD1InvB1, cvCholschurA1, cvDInvB_pn, cvDInvB_nrn, cvCholschurA);
            }

            R_chk_free(cvVz);
            R_chk_free(cvCholIplusXTildeVzXTildet);
            R_chk_free(tmp_n1n1r);
            if(failCode != 0){ goto cv_done; }

          }

          // Monte Carlo LOO-PD of the nEps models, sharing the pre-processing above
          for(e = 0; e < nEps; e++){

            epsilon = epsVec[e];

            // Fit on block-deleted data and obtain LOO-PD by Monte Carlo average
            int nBlockr = 0;
            double *Zb = NULL;

            for(sMC_CV = 0; sMC_CV < loopd_nMC; sMC_CV += nBlockMax){

              nBlock = std::min(nBlockMax, loopd_nMC - sMC_CV);

              // random variates of the block, in the order of a draw-by-draw loop
              for(bb = 0; bb < nBlock; bb++){


                if(family == family_poisson){
                  for(cv_i = 0; cv_i < nnk; cv_i++){
                    dtemp1 = cvY[cv_i] + epsilon;
                    dtemp2 = 1.0;
                    CV_eta[bb*nnk + cv_i] = rlogGamma(dtemp1);                              // log(Gamma(y + epsilon, 1)), underflow-safe
                  }
                }

                if(family == family_binomial){
                  for(cv_i = 0; cv_i < nnk; cv_i++){
                    dtemp1 = cvY[cv_i] + epsilon;
                    dtemp2 = cv_nBinom[cv_i];
                    dtemp2 += 2.0 * epsilon;
                    dtemp2 -= dtemp1;
                    CV_eta[bb*nnk + cv_i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
                  }
                }

                if(family == family_binary){
                  for(cv_i = 0; cv_i < nnk; cv_i++){
                    dtemp1 = cvY[cv_i] + epsilon;
                    dtemp2 = cv_nBinom[cv_i];
                    dtemp2 += 2.0 * epsilon;
                    dtemp2 -= dtemp1;
                    CV_eta[bb*nnk + cv_i] = rlogitBeta(dtemp1, dtemp2);                    // logit(Beta(y + epsilon, n - y + epsilon)), no rounding to 0 or 1
                  }
                }

                dtemp1 = 0.5 * nu_beta;
                dtemp2 = 1.0 / dtemp1;
                dtemp3 = rgamma(dtemp1, dtemp2);
                dtemp3 = 1.0 / dtemp3;
                dtemp3 = sqrt(dtemp3);
                for(j = 0; j < p; j++){
                  CV_beta[bb*p + j] = rnorm(0.0, dtemp3);
                }

                for(cv_i = 0; cv_i < nnk; cv_i++){
                  CV_xi[bb*nnk + cv_i] = rnorm(0.0, sigma_xi);
                }

                if(processType == "independent.shared"){
                  dtemp1 = 0.5 * nu_z;
                  dtemp2 = 1.0 / dtemp1;
                  dtemp3 = rgamma(dtemp1, dtemp2);
                  dtemp3 = 1.0 / dtemp3;
                  dtemp3 = sqrt(dtemp3);
                  for(k = 0; k < r; k++){
                    for(cv_i = 0; cv_i < nnk; cv_i++){
                      CV_z[bb*nnkr + k*nnk + cv_i] = rnorm(0.0, dtemp3);
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
                      CV_z[bb*nnkr + k*nnk + cv_i] = rnorm(0.0, dtemp3);
                    }
                  }
                }else if(processType == "multivariate"){

                  for(k = 0; k < r; k++){
                    for(cv_i = 0; cv_i < nnk; cv_i++){
                      tmp_n1r[k*nnk + cv_i] = rnorm(0.0, 1.0);
                    }
                  }
                  if(rInvWishart(r, nu_z + 2*r, chol_iwScale, samp_Sigma, tmp_rr) != 0){ failCode = 5; goto cv_done; }
                  F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){ failCode = 5; goto cv_done; }
                  mkLT(samp_Sigma, r);
                  F77_NAME(dgemm)(ntran, ytran, &nnk, &r, &r, &one, tmp_n1r, &nnk, samp_Sigma, &r, &zero, &CV_z[bb*nnkr], &nnk FCONE FCONE);

                }

                // variates of the z_tilde draws, scaled below once the projection gives the scales:
                // per process rgamma(0.5*(nu_z + nnk), .) and the nk standard normals of rnorm(0, sd) = sd*norm_rand()
                // (independent processes), or the Bartlett factor of the inverse-Wishart draw and nk x r standard
                // normals (multivariate)
                if(processType == "independent.shared" || processType == "independent"){
                  for(k = 0; k < r; k++){
                    dtemp1 = 0.5 * (nu_z + nnk);
                    dtemp2 = 1.0 / dtemp1;
                    CV_gam[bb*r + k] = rgamma(dtemp1, dtemp2);
                    for(cv_i = 0; cv_i < nk; cv_i++){
                      CV_nrm[bb*nkmaxr + k*nk + cv_i] = norm_rand();
                    }
                  }
                }else if(processType == "multivariate"){
                  rWishartBartlett(r, nu_z + nnk + 2*r, &CV_bart[bb*rr]);
                  for(k = 0; k < r; k++){
                    for(cv_i = 0; cv_i < nk; cv_i++){
                      CV_nrm[bb*nkmaxr + k*nk + cv_i] = norm_rand();                                       // = rnorm(0.0, 1.0)
                    }
                  }
                }

              }

              // projection step for the block
              projGLMvcbatch(nnk, p, r, nBlock, cvX, cvX_tilde, sigmaSq_xi, Lbeta, cvCholVz, processType,
                             CV_eta, CV_xi, CV_beta, CV_z, cvD1Inv, cvD1InvB1, cvCholschurA1,
                             cvDInvB_pn, cvDInvB_nrn, cvCholschurA, CV_tmpnr, CV_tmpn, CV_tmpp);

              // CV_z <- inv(Lz)*CV_z for every process block of every draw
              nBlockr = nBlock * r;
              if(processType == "independent.shared" || processType == "multivariate"){
                F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nBlockr, &one, cvCholVz, &nnk, CV_z, &nnk FCONE FCONE FCONE FCONE);
              }else if(processType == "independent"){
                for(k = 0; k < r; k++){
                  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nnk, &nBlock, &one, &cvCholVz[nnknnk * k], &nnk, &CV_z[nnk * k], &nnkr FCONE FCONE FCONE FCONE);
                }
              }

              // Prediction at held-out points for each draw of the block
              for(bb = 0; bb < nBlock; bb++){

                Zb = &CV_z[bb*nnkr];

                if(processType == "independent.shared"){

                  F77_NAME(dgemm)(ytran, ntran, &nk, &r, &nnk, &one, cvCz, &nnk, Zb, &nnk, &zero, z_tilde_mu, &nk FCONE FCONE);   // z_tilde_mu <- t(Cz)*inv(Vz)*v_z
                  for(k = 0; k < r; k++){
                    PCM_dist[0] = pow(F77_NAME(dnrm2)(&nnk, &Zb[nnk * k], &incOne), 2);

                    // sample z_tilde
                    dtemp2 = 1.0 / CV_gam[bb*r + k];
                    dtemp1 = (PCM_dist[0] + nu_z) / (nu_z + nnk);
                    dtemp3 = dtemp1 * dtemp2;
                    dtemp1 = sqrt(dtemp3);
                    for(cv_i = 0; cv_i < nk; cv_i++){
                      z_tilde[k*nk + cv_i] = dtemp1 * CV_nrm[bb*nkmaxr + k*nk + cv_i];                                 // = rnorm(0.0, dtemp1)
                    }
                  }
                  F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
                  F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, tmp_nkmaxr, &incOne);
                  F77_NAME(dcopy)(&nkr, tmp_nkmaxr, &incOne, z_tilde, &incOne);

                }else if(processType == "independent"){

                  for(k = 0; k < r; k++){
                    F77_NAME(dgemv)(ytran, &nnk, &nk, &one, &cvCz[nnk * nk * k], &nnk, &Zb[nnk * k], &incOne, &zero, &z_tilde_mu[nk * k], &incOne FCONE);
                    PCM_dist[0] = pow(F77_NAME(dnrm2)(&nnk, &Zb[nnk * k], &incOne), 2);

                    // sample z_tilde
                    dtemp2 = 1.0 / CV_gam[bb*r + k];
                    dtemp1 = (PCM_dist[0] + nu_z) / (nu_z + nnk);
                    dtemp3 = dtemp1 * dtemp2;
                    dtemp1 = sqrt(dtemp3);
                    for(cv_i = 0; cv_i < nk; cv_i++){
                      z_tilde[k*nk + cv_i] = dtemp1 * CV_nrm[bb*nkmaxr + k*nk + cv_i];                                 // = rnorm(0.0, dtemp1)
                    }
                    F77_NAME(dgemv)(ntran, &nk, &nk, &one, &z_tilde_cov[nknk * k], &nk, &z_tilde[nk * k], &incOne, &zero, &tmp_nkmaxr[nk * k], &incOne FCONE);
                  }
                  F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, tmp_nkmaxr, &incOne);
                  F77_NAME(dcopy)(&nkr, tmp_nkmaxr, &incOne, z_tilde, &incOne);

                }else if(processType == "multivariate"){

                  F77_NAME(dgemm)(ytran, ntran, &nk, &r, &nnk, &one, cvCz, &nnk, Zb, &nnk, &zero, z_tilde_mu, &nk FCONE FCONE);   // z_tilde_mu <- t(C)*inv(R)*v_z
                  F77_NAME(dgemm)(ytran, ntran, &r, &r, &nnk, &one, Zb, &nnk, Zb, &nnk, &zero, PCM_dist, &r FCONE FCONE);          // PCM_dist <- t(v_z)*inv(R)*v_z
                  F77_NAME(daxpy)(&rr, &one, iwScale, &incOne, PCM_dist, &incOne);                                                 // PCM_dist <- iwScale + t(v_z)*inv(R)*v_z
                  F77_NAME(dpotrf)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){ failCode = 5; goto cv_done; }
                  F77_NAME(dpotri)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){ failCode = 5; goto cv_done; }
                  F77_NAME(dpotrf)(lower, &r, PCM_dist, &r, &info FCONE); if(info != 0){ failCode = 5; goto cv_done; }
                  mkLT(PCM_dist, r);
                  if(invWishartFromBartlett(r, &CV_bart[bb*rr], PCM_dist, samp_Sigma, tmp_rr) != 0){ failCode = 5; goto cv_done; }  // = rInvWishart(r, nu_z + nnk + 2r, PCM_dist, ...)
                  F77_NAME(dpotrf)(lower, &r, samp_Sigma, &r, &info FCONE); if(info != 0){ failCode = 5; goto cv_done; }
                  mkLT(samp_Sigma, r);

                  for(k = 0; k < r; k++){
                    for(cv_i = 0; cv_i < nk; cv_i++){
                      z_tilde[k*nk + cv_i] = CV_nrm[bb*nkmaxr + k*nk + cv_i];
                    }
                  }
                  F77_NAME(dgemm)(ntran, ntran, &nk, &r, &nk, &one, z_tilde_cov, &nk, z_tilde, &nk, &zero, tmp_nkmaxr, &nk FCONE FCONE);
                  F77_NAME(dgemm)(ntran, ytran, &nk, &r, &r, &one, tmp_nkmaxr, &nk, samp_Sigma, &r, &zero, z_tilde, &nk FCONE FCONE);
                  F77_NAME(daxpy)(&nkr, &one, z_tilde_mu, &incOne, z_tilde, &incOne);

                }

                // Find canonical parameter (X*beta + X_tilde*z_tilde)
                lmulm_XTilde_VC(ntran, nk, r, 1, X_tilde_pred, z_tilde, tmp_nkmaxr);
                F77_NAME(dgemv)(ntran, &nk, &p, &one, X_pred, &nk, &CV_beta[bb*p], &incOne, &one, tmp_nkmaxr, &incOne FCONE);

                // Find CV-LOO-PD
                if(family == family_poisson){
                  for(cv_i = 0; cv_i < nk; cv_i++){
                    dtemp1 = exp(tmp_nkmaxr[cv_i]);
                    loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV + bb] = dpois(Y_pred[cv_i], dtemp1, 1);
                  }
                }

                if(family == family_binomial){
                  for(cv_i = 0; cv_i < nk; cv_i++){
                    dtemp1 = inverse_logit(tmp_nkmaxr[cv_i]);
                    loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV + bb] = dbinom(Y_pred[cv_i], nBinom_pred[cv_i], dtemp1, 1);
                  }
                }

                if(family == family_binary){
                  for(cv_i = 0; cv_i < nk; cv_i++){
                    dtemp1 = inverse_logit(tmp_nkmaxr[cv_i]);
                    loopd_val_MC_CV[cv_i*loopd_nMC + sMC_CV + bb] = dbinom(Y_pred[cv_i], 1.0, dtemp1, 1);
                  }
                }

              }

            }

            for(cv_i = 0; cv_i < nk; cv_i++){
              REAL(VECTOR_ELT(loopd_out_l, e))[start_index + cv_i] = logMeanExp(&loopd_val_MC_CV[cv_i*loopd_nMC], loopd_nMC);
            }

          }

        }

        cv_done:

        PutRNGstate();

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
        R_chk_free(cvCholVz);
        R_chk_free(cvCz);
        R_chk_free(z_tilde_cov);
        R_chk_free(z_tilde_mu);
        R_chk_free(z_tilde);
        R_chk_free(PCM_dist);
        R_chk_free(cvD1Inv);
        R_chk_free(cvD1InvB1);
        R_chk_free(cvCholschurA1);
        R_chk_free(cvDInvB_pn);
        R_chk_free(cvDInvB_nrn);
        R_chk_free(cvCholschurA);
        R_chk_free(del_PB);
        R_chk_free(del_QB);
        R_chk_free(del_QBK);
        R_chk_free(del_LP);
        R_chk_free(del_LQ);
        R_chk_free(del_WB);
        R_chk_free(del_Z);
        R_chk_free(del_H);
        R_chk_free(del_HK);
        R_chk_free(del_A2);
        R_chk_free(del_tmp_np);
        R_chk_free(del_w);
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
        R_chk_free(CV_eta);
        R_chk_free(CV_xi);
        R_chk_free(CV_tmpn);
        R_chk_free(CV_z);
        R_chk_free(CV_tmpnr);
        R_chk_free(CV_beta);
        R_chk_free(CV_tmpp);
        R_chk_free(CV_gam);
        R_chk_free(CV_nrm);
        R_chk_free(CV_bart);

      }

    }

    if(failCode != 0){
      R_chk_free(D1Inv); R_chk_free(D1InvB1); R_chk_free(DInvB_pn); R_chk_free(DInvB_nrn); R_chk_free(cholschurA);
      UNPROTECT(nProtect);
      glmLOOError(failCode);
    }

    // return object: a list of nEps fits, each a list of the posterior samples (and leave-one-out predictive densities)
    int nResultListObjs = loopd ? 4 : 3;
    SEXP fit_r;
    result_r = PROTECT(Rf_allocVector(VECSXP, nEps)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("z"));
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("xi"));
    if(loopd){
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("loopd"));
    }

    for(e = 0; e < nEps; e++){
      fit_r = Rf_allocVector(VECSXP, nResultListObjs);
      SET_VECTOR_ELT(result_r, e, fit_r);                                                // protected through result_r
      SET_VECTOR_ELT(fit_r, 0, VECTOR_ELT(samples_beta_l, e));                           // samples of beta
      SET_VECTOR_ELT(fit_r, 1, VECTOR_ELT(samples_z_l, e));                              // samples of z
      SET_VECTOR_ELT(fit_r, 2, VECTOR_ELT(samples_xi_l, e));                             // samples of xi
      if(loopd){
        SET_VECTOR_ELT(fit_r, 3, VECTOR_ELT(loopd_out_l, e));                            // leave-one-out predictive densities
      }
      Rf_namesgets(fit_r, resultName_r);
      SET_VECTOR_ELT(result_r, e, appendDiagnosticsRows(fit_r, diagPivot, diagMinCor, diagMaxCor, nDiag));
    }

    R_chk_free(D1Inv);
    R_chk_free(D1InvB1);
    R_chk_free(DInvB_pn);
    R_chk_free(DInvB_nrn);
    R_chk_free(cholschurA);

    UNPROTECT(nProtect);
    // return R_NilValue;

    return result_r;

}

extern "C" {

  // fit of a single candidate model
  SEXP stvcGLMexactLOO(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                       SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                       SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                       SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                       SEXP CV_K_r, SEXP loopd_nMC_r, SEXP cvUpdate_r, SEXP verbose_r){

    if(Rf_length(epsilon_r) != 1){
      Rf_error("c++ error: stvcGLMexactLOO expects a single boundary adjustment parameter.");
    }

    SEXP result_r = PROTECT(stvcGLMexactLOO_fit(Y_r, X_r, X_tilde_r, n_r, p_r, r_r, family_r, nBinom_r,
                                                sp_coords_r, time_coords_r, corfn_r,
                                                betaV_r, nu_beta_r, nu_z_r, sigmaSq_xi_r, iwScale_r,
                                                processType_r, phi_s_r, phi_t_r, epsilon_r,
                                                nSamples_r, loopd_r, loopd_method_r,
                                                CV_K_r, loopd_nMC_r, cvUpdate_r, verbose_r));
    UNPROTECT(1);

    return VECTOR_ELT(result_r, 0);

  } // end stvcGLMexactLOO

  // Fits of the candidate models sharing (phi_s, phi_t) for a vector of boundary adjustment parameters epsilon_r.
  // All the pre-processing depends on (phi_s, phi_t) only; it is computed once and shared by the fits (see
  // stvcGLMexactLOO_fit). Returns a list with one element per epsilon, each as returned by stvcGLMexactLOO.
  SEXP stvcGLMexactLOOgrid(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                           SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                           SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                           SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                           SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                           SEXP CV_K_r, SEXP loopd_nMC_r, SEXP cvUpdate_r, SEXP verbose_r){

    return stvcGLMexactLOO_fit(Y_r, X_r, X_tilde_r, n_r, p_r, r_r, family_r, nBinom_r,
                               sp_coords_r, time_coords_r, corfn_r,
                               betaV_r, nu_beta_r, nu_z_r, sigmaSq_xi_r, iwScale_r,
                               processType_r, phi_s_r, phi_t_r, epsilon_r,
                               nSamples_r, loopd_r, loopd_method_r,
                               CV_K_r, loopd_nMC_r, cvUpdate_r, verbose_r);

  } // end stvcGLMexactLOOgrid

}