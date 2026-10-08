#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
#include "psis.h"
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

  SEXP spLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                    SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                    SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                    SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
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

    // Set-up distance matrix and spatial correlation function
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

    // Leave-one-out predictive density details
    int loopd = INTEGER(loopd_r)[0];
    std::string loopd_method = CHAR(STRING_ELT(loopd_method_r, 0));
    const char *exact_str = "exact";
    const char *psis_str = "psis";

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
        Rprintf("\tbeta: Gaussian\n");
        Rprintf("\tmu:"); printVec(betaMu, p);
        Rprintf("\tcov:\n"); printMtrx(betaV, p, p);
        Rprintf("\n");
      }

      // prior on the measurement error variance sigma.sq; b = 0 corresponds
      // to the flat prior p(sigma.sq) proportional to 1/sigma.sq
      if(sigmaSqIGb == 0.0){
        Rprintf("\tsigma.sq: flat, proportional to 1/sigma.sq.\n\n");
      }else{
        Rprintf("\tsigma.sq: Inverse-Gamma\n\tshape = %.2f, scale = %.2f.\n\n",
                sigmaSqIGa, sigmaSqIGb);
      }

      Rprintf("Spatial process parameters:\n");

      if(corfn == "matern"){
        Rprintf("\tphi = %.2f, and, nu = %.2f.\n", phi, nu);
      }else{
        Rprintf("\tphi = %.2f.\n", phi);
      }
      Rprintf("Noise-to-spatial variance ratio = %.2f.\n\n", deltasq);

      Rprintf("Number of posterior samples = %i.\n\n", nSamples);

      if(loopd){
        Rprintf("LOO-PD calculation method = %s.\n", loopd_method.c_str());
      }

      Rprintf("----------------------------------------\n");

    }

    /*****************************************
     Priors on (beta, sigmaSqz) for sampling
     *****************************************/
    // The sampler below is written in terms of the spatial variance
    // sigmaSqz = sigmaSq/deltasq, where sigmaSq is the measurement error
    // variance. The prior sigmaSq ~ IG(a, b), beta | sigmaSq ~ N(muBeta, sigmaSq*Vbeta)
    // is equivalent to sigmaSqz ~ IG(a, b/deltasq), beta | sigmaSqz ~ N(muBeta, sigmaSqz*deltasq*Vbeta).
    // The flat prior p(sigmaSq) proportional to 1/sigmaSq (a = b = 0) is
    // invariant to this rescaling and so is the flat prior on beta.
    sigmaSqIGb = sigmaSqIGb / deltasq;                                                     // b = b/deltasq
    if(betaPrior == "normal"){
      F77_NAME(dscal)(&pp, &deltasq, betaV, &incOne);                                      // betaV = deltasq*Vbeta
    }

    /*****************************************
     Set-up posterior sample vector/matrices etc.
     *****************************************/
    double sigmaSqIGaPost = 0, sigmaSqIGbPost = 0;
    double sse = 0;
    double dtemp = 0;

    // const double deltasqInv = 1.0 / deltasq;
    const double delta = sqrt(deltasq);

    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);              // correlation matrix
    double *cholVy = (double *) R_alloc(nn, sizeof(double)); zeros(cholVy, nn);      // allocate memory for n x n matrix
    double *thetasp = (double *) R_alloc(2, sizeof(double));                         // spatial process parameters

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);          // allocate memory for n x 1 vector

    double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);                  // allocate memory for p x 1 vector
    double *VbetaInvMuBeta = (double *) R_alloc(p, sizeof(double)); zeros(VbetaInvMuBeta, p);  // allocate memory for p x 1 vector

    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);  // allocate VbetaInv
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);      // allocate memory for p x p matrix
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp2, pp);    // allocate memory for p x p matrix

    //construct covariance matrix (full)
    thetasp[0] = phi;
    thetasp[1] = nu;
    spCorFull(coordsD, n, thetasp, corfn, Vz);

    // construct marginal covariance matrix (Vz+deltasq*I)
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVy, &incOne);
    for(i = 0; i < n; i++){
      cholVy[i*n + i] += deltasq;
    }

    // chol(Vy)
    F77_NAME(dpotrf)(lower, &n, cholVy, &n, &info FCONE); if(info != 0){perror("c++ error: Vy dpotrf failed\n");}

    // find cholinv(Vy)*Y
    F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);                                         // tmp_n = Y
    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVy, &n, tmp_n, &incOne FCONE FCONE FCONE);  // tmp_n = cholinv(Vy)*Y

    if(betaPrior == "normal"){
      // find VbetaInvmuBeta
      F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                     // VbetaInv = Vbeta
      F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");} // VbetaInv = chol(Vbeta)
      F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");} // VbetaInv = chol2inv(Vbeta)
      F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betaMu, &incOne, &zero, VbetaInvMuBeta, &incOne FCONE);       // VbetaInvMuBeta = VbetaInv*muBeta
    }else{
      // flat prior on beta: VbetaInv = 0 and VbetaInv*muBeta = 0 (already
      // zero-initialized); since p(beta) does not carry the factor
      // (sigmaSq)^(-p/2) of the conjugate prior, the posterior inverse-gamma
      // shape is aIG + (n - p)/2
      sigmaSqIGa -= 0.5 * p;
    }

    //  find XtVyInvY
    double *tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np, np);                            // allocate temporary memory for n x p matrix
    F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);                                                          // tmp_np = X
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, cholVy, &n, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = cholinv(Vy)*X
    F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);             // tmp_p1 = t(X)*VyInv*Y

    // find betahat = inv(XtVyInvX + VbetaInv)(XtVyInvY + VbetaInvmuBeta)
    F77_NAME(daxpy)(&p, &one, VbetaInvMuBeta, &incOne, tmp_p1, &incOne);                                        // tmp_p1 = XtVyInvY + VbetaInvmuBeta
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);     // tmp_pp = t(X)*VyInv*X

    F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, tmp_pp, &incOne);                                             // tmp_pp = t(X)*VyInv*X + VbetaInv
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}  // tmp_pp = chol(XtVyInvX + VbetaInv)
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                    // tmp_p1 = cholinv(XtVyInvX + VbetaInv)*tmp_p1

    // find sse = t(Y-X*betahat)*VyInv*(Y-X*betahat) + t(betahat-muBeta)*VbetaInv*(betahat-muBeta);
    // equal to t(Y)*VyInv*Y + t(muBeta)*VbetaInv*muBeta - t(m)*M*m, but each term is a sum of
    // squares, hence non-negative and free of catastrophic cancellation
    double *betahat = (double *) R_chk_calloc(p, sizeof(double)); zeros(betahat, p);                           // allocate temporary memory for p x 1 vector
    F77_NAME(dcopy)(&p, tmp_p1, &incOne, betahat, &incOne);                                                     // betahat = cholinv(XtVyInvX + VbetaInv)*tmp_p1
    F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, betahat, &incOne FCONE FCONE FCONE);                   // betahat = inv(XtVyInvX + VbetaInv)*(XtVyInvY + VbetaInvmuBeta)
    F77_NAME(dgemv)(ntran, &n, &p, &negOne, tmp_np, &n, betahat, &incOne, &one, tmp_n, &incOne FCONE);          // tmp_n = cholinv(Vy)*(Y-X*betahat)
    sse = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne), 2);                                                          // sse = t(Y-X*betahat)*VyInv*(Y-X*betahat)

    // deallocate tmp_np
    R_chk_free(tmp_np);

    if(betaPrior == "normal"){
      F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, betahat, &incOne);                                          // betahat = betahat - muBeta
      F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betahat, &incOne, &zero, tmp_n, &incOne FCONE);            // tmp_n[1:p] = VbetaInv*(betahat - muBeta)
      sse += F77_CALL(ddot)(&p, betahat, &incOne, tmp_n, &incOne);                                              // sse = sse + t(betahat-muBeta)*VbetaInv*(betahat-muBeta)
    }

    // deallocate betahat
    R_chk_free(betahat);

    if(!(sigmaSqIGb + 0.5 * sse > 0.0)){
      Rf_error("c++ error: posterior scale of sigmaSq is not positive (sse = %g).\n", sse);
    }

    // set-up for sampling spatial random effects
    double *tmp_nn2 = (double *) R_chk_calloc(nn, sizeof(double)); zeros(tmp_nn2, nn);                           // calloc n x n matrix
    F77_NAME(dcopy)(&nn, Vz, &incOne, tmp_nn2, &incOne);                                                         // tmp_nn2 = Vz
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &n, &one, cholVy, &n, tmp_nn2, &n FCONE FCONE FCONE FCONE);  // tmp_nn2 = cholinv(Vy)*Vz
    F77_NAME(dtrsm)(lside, lower, ytran, nUnit, &n, &n, &one, cholVy, &n, tmp_nn2, &n FCONE FCONE FCONE FCONE);  // tmp_nn2 = inv(Vy)*Vz
    F77_NAME(dpotrf)(lower, &n, tmp_nn2, &n, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}  // tmp_nn2 = chol(inv(Vy)*Vz)
    mkLT(tmp_nn2, n);                                                                                            // make cholDinv lower-triangular

    // posterior parameters of sigmaSq
    sigmaSqIGaPost += sigmaSqIGa;
    sigmaSqIGaPost += 0.5 * n;

    sigmaSqIGbPost += sigmaSqIGb;
    sigmaSqIGbPost += 0.5 * sse;

    // posterior samples of sigma-sq and beta
    SEXP samples_sigmaSq_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_sigmaSqz_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

    // sample storage at s-th iteration temporary allocation
    double sigmaSqz = 0;
    double *beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta, p);
    double *z = (double *) R_chk_calloc(n, sizeof(double)); zeros(z, n);

    GetRNGstate();

    for(s = 0; s < nSamples; s++){
      // sample sigmaSqz (spatial variance) from its marginal posterior
      dtemp = 1.0 / sigmaSqIGbPost;
      dtemp = rgamma(sigmaSqIGaPost, dtemp);
      sigmaSqz = 1.0 / dtemp;
      REAL(samples_sigmaSqz_r)[s] = sigmaSqz;
      REAL(samples_sigmaSq_r)[s] = deltasq * sigmaSqz;                                       // sigmaSq = deltasq*sigmaSqz

      // sample fixed effects by composition sampling
      dtemp = sqrt(sigmaSqz);
      for(j = 0; j < p; j++){
        beta[j] = rnorm(tmp_p1[j], dtemp);                                                   // beta ~ N(tmp_p1, sigmaSq*I)
      }
      F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, beta, &incOne FCONE FCONE FCONE); // beta = t(cholinv(tmp_pp))*beta

      dtemp = dtemp * delta;                                                                 // dtemp = sqrt(deltasq*sigmaSq)
      // sample spatial effects by composition sampling
      for(i = 0; i < n; i++){
        tmp_n[i] = rnorm(0.0, dtemp);                                                               // tmp_n ~ N(0, sigmaSq*I)
      }
      F77_NAME(dcopy)(&n, Y, &incOne, z, &incOne);                                                  // z = Y
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &one, z, &incOne FCONE);        // z = Y-X*beta
      F77_NAME(dgemv)(ytran, &n, &n, &one, tmp_nn2, &n, z, &incOne, &one, tmp_n, &incOne FCONE);    // tmp_n = tmp_n + t(chol(tmp_nn2))*(Y-X*beta)/deltasq
      F77_NAME(dgemv)(ntran, &n, &n, &one, tmp_nn2, &n, tmp_n, &incOne, &zero, z, &incOne FCONE);   // z = chol(tmp_nn2)*tmp_n

      // copy samples into SEXP return object
      F77_NAME(dcopy)(&p, &beta[0], &incOne, &REAL(samples_beta_r)[s*p], &incOne);
      F77_NAME(dcopy)(&n, &z[0], &incOne, &REAL(samples_z_r)[s*n], &incOne);

    }

    PutRNGstate();

    // Free stuff
    R_chk_free(tmp_nn2);
    R_chk_free(beta);
    R_chk_free(z);

    // make return object
    SEXP result_r, resultName_r;

    // If loopd is TRUE, set-up Leave-one-out predictive density calculation
    if(loopd){

      int n1 = n - 1;
      int n1n1 = n1 * n1;
      int n1p = n1 * p;

      SEXP loopd_out_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;
      SEXP loopd_k_r = R_NilValue;                                    // Pareto k diagnostics (PSIS only)

      // Exact leave-one-out predictive densities calculation
      if(loopd_method == exact_str){

        const int marginal_model = 1;        // SWITCH = 1: Use marginal model to find LOO-PD; 0: Use unmarginalized augmented model

        if(marginal_model){

          double *looX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looX, n1p);
          double *X_tilde = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_tilde, p);
          double *looY = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looY, n1);
          double *looJ = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looJ, n1);
          double *looCholVy = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVy, n1n1);
          double *looB = (double *) R_chk_calloc(pp, sizeof(double)); zeros(looB, pp);
          double *looBb = (double *) R_chk_calloc(p, sizeof(double)); zeros(looBb, p);
          double *looH = (double *) R_chk_calloc(p, sizeof(double)); zeros(looH, p);
          double *tmp_n11 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n11, n1);

          double location = 0.0, scale = 0.0, loo_sse = 0.0, a_star = 0.0, b_star = 0.0;

          a_star = sigmaSqIGa;
          a_star += 0.5 * n1;

          int loo_index = 0;

          for(loo_index = 0; loo_index < n; loo_index++){

            copyMatrixDelRow(X, n, p, looX, loo_index);                   // Row-deleted X
            copyMatrixRowToVec(X, n, p, X_tilde, loo_index);              // Copy left out X into Xtilde
            copyVecExcludingOne(Y, looY, n, loo_index);                   // Leave-one-out Y
            cholRowDelUpdate(n, cholVy, loo_index, looCholVy, tmp_n11);   // Row-deletion Cholesky update of Vy
            copyVecExcludingOne(&Vz[loo_index*n], looJ, n, loo_index);    // h2 = Vz[-i,i]

            F77_NAME(dtrsv)(lower, ntran, nUnit, &n1, looCholVy, &n1, looY, &incOne FCONE FCONE FCONE); // looY = cholinv(looVy)*Y[-i]
            F77_NAME(dtrsv)(lower, ntran, nUnit, &n1, looCholVy, &n1, looJ, &incOne FCONE FCONE FCONE); // looJ = cholinv(looVy)*J

            location = F77_CALL(ddot)(&n1, looJ, &incOne, looY, &incOne);                               // location = t(J)*inv(Vy)*Y

            F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n1, &p, &one, looCholVy, &n1, looX, &n1 FCONE FCONE FCONE FCONE);  // looX = cholinv(looVy)*looX
            F77_NAME(dgemm)(ytran, ntran, &p, &p, &n1, &one, looX, &n1, looX, &n1, &zero, looB, &p FCONE FCONE);            // looB = t(X)*VyInv*X
            F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, looB, &incOne);                                                   // looB = t(X)*VyInv*X + VbetaInv
            F77_NAME(dpotrf)(lower, &p, looB, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}        // looB = chol(t(X)*VyInv*X + VbetaInv)
            F77_NAME(dgemv)(ytran, &n1, &p, &one, looX, &n1, looY, &incOne, &zero, looBb, &incOne FCONE);                   // looBb = t(X)*VyInv*Y
            F77_NAME(daxpy)(&p, &one, VbetaInvMuBeta, &incOne, looBb, &incOne);                                             // looBb = XtVyInvY + VbetaInvmuBeta
            F77_NAME(dtrsv)(lower, ntran, nUnit, &p, looB, &p, looBb, &incOne FCONE FCONE FCONE);                           // looBb = cholinv(looB)*b
            F77_NAME(dtrsv)(lower, ytran, nUnit, &p, looB, &p, looBb, &incOne FCONE FCONE FCONE);                           // looBb = inv(looB)*b = Bb

            // loo_sse in residual form (see sse above): no cancellation, always non-negative
            F77_NAME(dgemv)(ntran, &n1, &p, &negOne, looX, &n1, looBb, &incOne, &one, looY, &incOne FCONE);                 // looY = cholinv(looVy)*(Y[-i]-X[-i,]*Bb)
            loo_sse = pow(F77_NAME(dnrm2)(&n1, looY, &incOne), 2);                                                          // loo_sse = t(Y[-i]-X[-i,]*Bb)*looVyInv*(Y[-i]-X[-i,]*Bb)
            if(betaPrior == "normal"){
              F77_NAME(dcopy)(&p, looBb, &incOne, looH, &incOne);                                                           // looH = Bb
              F77_NAME(daxpy)(&p, &negOne, betaMu, &incOne, looH, &incOne);                                                 // looH = Bb - muBeta
              F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, looH, &incOne, &zero, tmp_n11, &incOne FCONE);                 // tmp_n11[1:p] = VbetaInv*(Bb - muBeta)
              loo_sse += F77_CALL(ddot)(&p, looH, &incOne, tmp_n11, &incOne);                                               // loo_sse = loo_sse + t(Bb-muBeta)*VbetaInv*(Bb-muBeta)
            }

            F77_NAME(dcopy)(&p, X_tilde, &incOne, looH, &incOne);                                                           // looH = X_tilde
            F77_NAME(dgemv)(ytran, &n1, &p, &negOne, looX, &n1, looJ, &incOne, &one, looH, &incOne FCONE);                  // looH = X_tilde - t(J)*VyInv*X

            location += F77_CALL(ddot)(&p, looH, &incOne, looBb, &incOne);                                                  // location = t(J)*inv(Vy)*Y + H*Bb

            scale = Vz[loo_index * n + loo_index] + deltasq;                                                                // scale = V(y_tilde)
            scale -= F77_CALL(ddot)(&n1, looJ, &incOne, looJ, &incOne);                                                     // scale = V(y_tilde) - t(J)*inv(Vy)*J

            F77_NAME(dtrsv)(lower, ntran, nUnit, &p, looB, &p, looH, &incOne FCONE FCONE FCONE);                            // looH = cholinv(looB)*looH

            scale += F77_CALL(ddot)(&p, looH, &incOne, looH, &incOne);                                                      // scale = V(y_tilde) - t(J)*inv(Vy)*J + H*B*t(H)

            b_star = sigmaSqIGb + 0.5 * loo_sse;
            scale = (b_star / a_star) * scale;
            scale = sqrt(scale);

            dtemp = Y[loo_index] - location;
            dtemp = dtemp / scale;
            dtemp = Rf_dt(dtemp, 2 * a_star, 1);
            dtemp -= log(scale);

            REAL(loopd_out_r)[loo_index] = dtemp;

          }

          R_chk_free(looX);
          R_chk_free(X_tilde);
          R_chk_free(looY);
          R_chk_free(looJ);
          R_chk_free(looCholVy);
          R_chk_free(looB);
          R_chk_free(looBb);
          R_chk_free(looH);
          R_chk_free(tmp_n11);

          // End LOO-PD calculation using marginal model

        }else{

          double *looX = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(looX, n1p);
          double *looVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looVz, n1n1);
          double *looCholVy = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVy, n1n1);
          double *cholVz = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholVz, nn);
          double *looCholVz = (double *) R_chk_calloc(n1n1, sizeof(double)); zeros(looCholVz, n1n1);
          double *h1 = (double *) R_chk_calloc(p, sizeof(double)); zeros(h1, p);
          double *h2 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(h2, n1);
          double *tmp_n11 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n11, n1);
          double *tmp_n12 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(tmp_n12, n1);
          double *tmp_n1p1 = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(tmp_n1p1, n1p);
          double *tmp_n1p2 = (double *) R_chk_calloc(n1p, sizeof(double)); zeros(tmp_n1p2, n1p);
          double *out_p = (double *) R_chk_calloc(p, sizeof(double)); zeros(out_p, p);
          double *out_n1 = (double *) R_chk_calloc(n1, sizeof(double)); zeros(out_n1, n1);

          const double deltasqInv = 1.0 / deltasq;

          int loo_index = 0;
          double location = 0.0, scale = 0.0, a_star = 0.0, b_star = 0.0;

          a_star = sigmaSqIGa;
          a_star += 0.5 * n1;

          F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
          F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}

          for(loo_index = 0; loo_index < n; loo_index++){

            copyMatrixDelRow(X, n, p, looX, loo_index);                   // Row-deleted X
            copyMatrixDelRowCol(Vz, n, n, looVz, loo_index, loo_index);   // Row-column deleted Vz
            cholRowDelUpdate(n, cholVy, loo_index, looCholVy, tmp_n11);   // Row-deletion CHOL update Vy
            cholRowDelUpdate(n, cholVz, loo_index, looCholVz, tmp_n11);   // Row-deletion CHOL update Vz
            copyMatrixRowToVec(X, n, p, h1, loo_index);                   // h1 = X[i,1:p]
            copyVecExcludingOne(&Vz[loo_index*n], h2, n, loo_index);      // h2 = Vz[-i,i]

            inversionLM(looX, n1, p, deltasq, VbetaInv, looVz, looCholVy, h1, h2,
                        tmp_n11, tmp_n12, tmp_p1, tmp_pp, tmp_n1p1,
                        out_p, out_n1, 1);                                // call inversionLM with LOO = TRUE

            F77_NAME(dtrsv)(lower, ntran, nUnit, &n1, looCholVz, &n1, h2, &incOne FCONE FCONE FCONE);
            F77_NAME(dtrsv)(lower, ytran, nUnit, &n1, looCholVz, &n1, h2, &incOne FCONE FCONE FCONE);
            scale = F77_CALL(ddot)(&p, out_p, &incOne, h1, &incOne);
            scale += F77_CALL(ddot)(&n1, out_n1, &incOne, h2, &incOne);

            copyVecExcludingOne(Y, h2, n, loo_index);                      // h2 = Y[-i]
            F77_NAME(dscal)(&n1, &deltasqInv, h2, &incOne);                // h2 = Y[-i]/deltasq
            location = F77_CALL(ddot)(&n1, out_n1, &incOne, h2, &incOne);  // loc = t(h2)*Mstar*Y[-i]/deltasq

            F77_NAME(dgemv)(ytran, &n1, &p, &one, looX, &n1, h2, &incOne, &zero, h1, &incOne FCONE); // h1 = t(X)Y[-i]/deltasq
            F77_NAME(daxpy)(&p, &one, VbetaInvMuBeta, &incOne, h1, &incOne);                         // h1 = t(X)Y[-i]/deltasq + VbetaInvMuBeta
            location += F77_CALL(ddot)(&p, out_p, &incOne, h1, &incOne);                             // loc = t(h)*Mstar*(gamma_hat)

            b_star = pow(F77_NAME(dnrm2)(&n1, h2, &incOne), 2);           // b_star = t(Y[-i])*Y[-i]/(deltasq)^2
            b_star *= deltasq;                                            // b_star = t(Y[-i])*Y[-i]/deltasq

            inversionLM(looX, n1, p, deltasq, VbetaInv, looVz, looCholVy, h1, h2,
                        tmp_n11, tmp_n12, tmp_p1, tmp_pp, tmp_n1p1,
                        out_p, out_n1, 0);                                // call inversionLM with LOO = FALSE
            b_star -= F77_CALL(ddot)(&p, out_p, &incOne, h1, &incOne);
            b_star -= F77_CALL(ddot)(&n1, out_n1, &incOne, h2, &incOne);  // b_star = sse_star
            b_star *= 0.5;
            b_star += sigmaSqIGb;

            scale += deltasq;
            scale = (b_star / a_star) * scale;
            scale = sqrt(scale);

            dtemp = Y[loo_index] - location;
            dtemp = dtemp / scale;
            dtemp = Rf_dt(dtemp, 2*a_star, 1);
            dtemp -= log(scale);

            REAL(loopd_out_r)[loo_index] = dtemp;

          }

          R_chk_free(looX);
          R_chk_free(looVz);
          R_chk_free(looCholVy);
          R_chk_free(cholVz);
          R_chk_free(looCholVz);
          R_chk_free(h1);
          R_chk_free(h2);
          R_chk_free(tmp_n11);
          R_chk_free(tmp_n12);
          R_chk_free(tmp_n1p1);
          R_chk_free(tmp_n1p2);
          R_chk_free(out_p);
          R_chk_free(out_n1);


        }


      }

      if(loopd_method == psis_str){

        int loo_index = 0, s = 0;
        double theta_i = 0.0;

        // PSIS workspace, O(nSamples), allocated once for all observations
        int psis_L = psis_tail_length(nSamples);
        int psis_M = psis_gpd_grid_length(psis_L);

        double *X_i = (double *) R_chk_calloc(p, sizeof(double)); zeros(X_i, p);
        double *ll_i = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(ll_i, nSamples);
        double *lw_i = (double *) R_chk_calloc(nSamples, sizeof(double)); zeros(lw_i, nSamples);
        int *idx_i = (int *) R_chk_calloc(nSamples, sizeof(int)); zeros(idx_i, nSamples);
        double *xtail_i = (double *) R_chk_calloc(psis_L, sizeof(double)); zeros(xtail_i, psis_L);
        double *theta_gpd = (double *) R_chk_calloc(psis_M, sizeof(double)); zeros(theta_gpd, psis_M);
        double *ltheta_gpd = (double *) R_chk_calloc(psis_M, sizeof(double)); zeros(ltheta_gpd, psis_M);

        loopd_k_r = PROTECT(Rf_allocVector(REALSXP, n)); nProtect++;

        double *pointer_beta = REAL(samples_beta_r);
        double *pointer_z = REAL(samples_z_r);
        double *pointer_sigmaSq = REAL(samples_sigmaSq_r);

        for(loo_index = 0; loo_index < n; loo_index++){

          copyMatrixRowToVec(X, n, p, X_i, loo_index);                                          // X_i = X[i,1:p]

          // log-likelihood of the i-th observation at each posterior draw
          for(s = 0; s < nSamples; s++){
            theta_i = F77_CALL(ddot)(&p, X_i, &incOne, &pointer_beta[s*p], &incOne);            // theta_i = X_i * beta_s
            theta_i += pointer_z[n*s + loo_index];                                              // theta_i = X_i*beta_s + zi_s
            ll_i[s] = Rf_dnorm4(Y[loo_index], theta_i, sqrt(pointer_sigmaSq[s]), 1);           // sigmaSq_s is the measurement error variance
          }

          psis_loo(ll_i, nSamples, psis_L, lw_i, idx_i, xtail_i, theta_gpd, ltheta_gpd,
                   &REAL(loopd_out_r)[loo_index], &REAL(loopd_k_r)[loo_index]);

        }

        R_chk_free(X_i);
        R_chk_free(ll_i);
        R_chk_free(lw_i);
        R_chk_free(idx_i);
        R_chk_free(xtail_i);
        R_chk_free(theta_gpd);
        R_chk_free(ltheta_gpd);

      }

      // make return object for posterior samples and leave-one-out predictive densities
      int nResultListObjs = (loopd_method == psis_str) ? 6 : 5;

      result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
      resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

      // samples of beta
      SET_VECTOR_ELT(result_r, 0, samples_beta_r);
      SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

      // samples of sigma-sq
      SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
      SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

      // samples of sigma-sq-z
      SET_VECTOR_ELT(result_r, 2, samples_sigmaSqz_r);
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("sigmaSq.z"));

      // samples of z
      SET_VECTOR_ELT(result_r, 3, samples_z_r);
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z"));

      // leave-one-out predictive densities
      SET_VECTOR_ELT(result_r, 4, loopd_out_r);
      SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("loopd"));

      // Pareto k diagnostics of PSIS
      if(loopd_method == psis_str){
        SET_VECTOR_ELT(result_r, 5, loopd_k_r);
        SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("loopd.pareto_k"));
      }

      Rf_namesgets(result_r, resultName_r);

    }else{

      // make return object for posterior samples of sigma-sq, beta and z
      int nResultListObjs = 4;

      result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
      resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

      // samples of beta
      SET_VECTOR_ELT(result_r, 0, samples_beta_r);
      SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

      // samples of sigma-sq
      SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
      SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

      // samples of sigma-sq-z
      SET_VECTOR_ELT(result_r, 2, samples_sigmaSqz_r);
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("sigmaSq.z"));

      // samples of z
      SET_VECTOR_ELT(result_r, 3, samples_z_r);
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("z"));

      Rf_namesgets(result_r, resultName_r);

    }

    // SEXP result_r = PROTECT(Rf_allocMatrix(REALSXP, nSamples, p)); nProtect++;
    // SEXP result_r1 = PROTECT(Rf_allocVector(REALSXP, 1)); nProtect++;
    // SEXP tmp_nn_r = PROTECT(Rf_allocMatrix(REALSXP, n, n)); nProtect++;

    // for (i = 0; i < n; i++) {
    //   for (j = 0; j < n; j++) {
    //     REAL(tmp_nn_r)[i*n + j] = tmp_nn2[i*n + j];
    //   }
    // }

    // for(i = 0; i < n; i++){
    //   REAL(tmp_n_r)[i] = tmp_n[i];
    // }

    // REAL(result_r1)[0] = sigmaSqIGaPost;
    // REAL(result_r1)[1] = sigmaSqIGbPost;

    // REAL(result_r1)[0] = sse;

    UNPROTECT(nProtect);

    return result_r;

  }

}
