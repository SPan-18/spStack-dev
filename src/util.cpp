#define USE_FC_LEN_T
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

void mysolve(double *A, double *b, int n){

  int info = 0;
  char const *upper = "U";
  char const *trans = "T";
  char const *ntrans = "N";
  char const *nunit = "N";
  int incx = 1;     // Increment for x

  F77_NAME(dpotrf)(upper, &n, A, &n, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
  F77_NAME(dtrsv)(upper, trans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);
  F77_NAME(dtrsv)(upper, ntrans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);

}
