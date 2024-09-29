#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP R_cholRankOneUpdate(SEXP L_r, SEXP n_r, SEXP v_r, SEXP alpha_r, SEXP beta_r, SEXP lower_r){

    double *L = REAL(L_r);
    double *v = REAL(v_r);
    int n = INTEGER(n_r)[0];
    double alpha = REAL(alpha_r)[0];
    double beta = REAL(beta_r)[0];
    int lower = INTEGER(lower_r)[0];

    int nn = n * n;
    SEXP L1 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    double *L1_pointer = REAL(L1); zeros(L1_pointer, nn);

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);

    if(lower){
      cholRankOneUpdate(n, L, alpha, beta, v, L1_pointer, tmp_n);
    }else{
      upperTri_lowerTri(L, n);
      mkLT(L, n);
      cholRankOneUpdate(n, L, alpha, beta, v, L1_pointer, tmp_n);
    }

    UNPROTECT(1);
    return L1;

  }

  SEXP R_cholRowDelUpdate(SEXP L_r, SEXP n_r, SEXP row_r, SEXP lower_r){

    double *L = REAL(L_r);
    int n = INTEGER(n_r)[0];
    int row_del = INTEGER(row_r)[0];
    int lower = INTEGER(lower_r)[0];

    int nMinusOne = n - 1;
    int n1n1 = nMinusOne * nMinusOne;
    SEXP L1 = PROTECT(Rf_allocMatrix(REALSXP, nMinusOne, nMinusOne));
    double *L1_pointer = REAL(L1); zeros(L1_pointer, n1n1);

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);

    // match index value with R indexing system
    row_del -= 1;

    if(lower){
      cholRowDelUpdate(n, L, row_del, L1_pointer, tmp_n);
    }else{
      upperTri_lowerTri(L, n);
      mkLT(L, n);
      cholRowDelUpdate(n, L, row_del, L1_pointer, tmp_n);
    }

    UNPROTECT(1);
    return L1;

  }

  SEXP R_cholRowBlockDelUpdate(SEXP L_r, SEXP n_r, SEXP start_r, SEXP end_r, SEXP lower_r){

    double *L = REAL(L_r);
    int n = INTEGER(n_r)[0];
    int del_start = INTEGER(start_r)[0];
    int del_end = INTEGER(end_r)[0];
    int lower = INTEGER(lower_r)[0];

    int nk = n - del_end + del_start - 1;
    int nknk = nk * nk;
    SEXP L1 = PROTECT(Rf_allocMatrix(REALSXP, nk, nk));
    double *L1_pointer = REAL(L1); zeros(L1_pointer, nknk);

    // match index value with R indexing system
    del_start -= 1;
    del_end -= 1;

    double *tmp_nk = (double *) R_alloc(nk, sizeof(double)); zeros(tmp_nk, nk);
    double *tmp_nknk = (double *) R_alloc(nknk, sizeof(double)); zeros(tmp_nknk, nknk);

    if(lower){
      cholBlockDelUpdate(n, L, del_start, del_end, L1_pointer, tmp_nknk, tmp_nk);
    }else{
      upperTri_lowerTri(L, n);
      mkLT(L, n);
      cholBlockDelUpdate(n, L, del_start, del_end, L1_pointer, tmp_nknk, tmp_nk);
    }

    UNPROTECT(1);
    return L1;

  }
}
