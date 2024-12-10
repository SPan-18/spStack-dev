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
                    SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r,
                    SEXP sharedProcess_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                    SEXP nSamples_r, SEXP verbose_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
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
    int nr = n * r;
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

    // spatial-temporal process parameters: create spatial-temporal covariance matrices
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    if(corfn == "gneiting-decay"){
        if(sharedProcess_r){
            double phi_s = REAL(phi_s_r)[0];
            double phi_t = REAL(phi_t_r)[0];
            thetaspt[0] = phi_s;
            thetaspt[1] = phi_t;
            double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
            sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, Vz);
            printMtrx(Vz, n, n);
        }else{
            double *phi_s = REAL(phi_s_r);
            double *phi_t = REAL(phi_t_r);
        }
    }

    return(R_NilValue);

    } // end stvcGLMexact
}