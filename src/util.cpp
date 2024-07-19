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

void mysolveUT(double *A, double *b, int n){

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

void mysolveLT(double *A, double *b, int n){

  int info = 0;
  char const *lower = "L";
  char const *trans = "T";
  char const *ntrans = "N";
  char const *nunit = "N";
  int incx = 1;     // Increment for x

  F77_NAME(dpotrf)(lower, &n, A, &n, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
  F77_NAME(dtrsv)(lower, ntrans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, trans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);

}

void printMtrx(double *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("%.3f\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }
}

void printVec(double *m, int n){

  Rprintf("\t");
  for(int j = 0; j < n; j++){
    Rprintf("%.3f\t", m[j]);
  }
  Rprintf("\n");
}

void printVec(int *m, int n){

  Rprintf(" ");
  for(int j = 0; j < n; j++){
    Rprintf("%i ", m[j]);
  }
  Rprintf("\n");
}

void spCorLT(double *D, int n, double *theta, std::string &corfn, double *C){
  int i,j;

  if(corfn == "exponential"){

    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){
        C[i*n + j] = theta[0] * exp(-1.0 * theta[1] * D[i*n + j]);
      }
    }

  }else if(corfn == "matern"){

    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){
        if(D[i*n + j] * theta[0] > 0.0){
          C[i*n + j] = pow(D[i*n + j] * theta[0], theta[1]) / (pow(2, theta[1] - 1) * gammafn(theta[1])) * bessel_k(D[i*n + j] * theta[0], theta[1], 1.0);
        }else{
          C[i*n + j] = 1.0;
        }
      }
    }

  }else{
    error("c++ error: corfn is not correctly specified");
  }
}

void zeros(double *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0.0;
}

void zeros(int *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0;
}
