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
                       SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r,
                       SEXP sharedProcess_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                       SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                       SEXP CV_K_r, SEXP loopd_nMC_r, SEXP verbose_r){

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

    return R_NilValue;

  }

}