#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP mysolveC(SEXP A_r, SEXP b_r, SEXP n_r){

    double *A = REAL(A_r);
    double *b = REAL(b_r);
    int n = INTEGER(n_r)[0];

    SEXP result = PROTECT(allocVector(REALSXP, n));

    mysolve(A, b, n);

    std::copy(b, b + n, REAL(result));

    UNPROTECT(1);
    return result;
  }
}
