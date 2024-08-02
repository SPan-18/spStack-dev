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

// Cholesky factor rank-1 update; chol(alpha*LLt + beta*vvt)
// void cholRankOneUpdate(int n, double *L, double alpha, double beta, double *v){
//
//
//
// }

// Copy matrix from C to SEXP
void copyMatrixSEXP(double *matrixC, int dim1, int dim2, double *pointerSEXP){

  int i, j;

  for(i = 0; i < dim2; i++){
    for(j = 0; j < dim1; j++){
      pointerSEXP[i*dim1 + j] = matrixC[i*dim1 + j];
    }
  }

}

// Copy vector from C to SEXP
void copyVectorSEXP(double *vectorC, int dim, double *pointerSEXP){

  int i;

  for(i = 0; i < dim; i++){
    pointerSEXP[i] = vectorC[i];
  }

}

// Copy a submatrix of A into a submatrix of B
void copySubmat(double *A, int nRowA, int nColA, double *B, int nRowB, int nColB,
                int startRowA, int startColA, int startRowB, int startColB,
                int nRowCopy, int nColCopy){

  if(startRowA + nRowCopy > nRowA || startColA + nColCopy > nColA){
    perror("Indices of rows/columns to copy exceeds dimensions of source matrix.");
  }

  if(startRowB + nRowCopy > nRowB || startColB + nColCopy > nColB){
    perror("Indices rows/columns to copy exceeds dimensions of destination matrix.");
  }

  int col, row;

  for(col = 0; col < nColCopy; col++){
    for(row= 0; row < nRowCopy; row++){
      B[(startColB + col)*nRowB + (startRowB + row)] = A[(startColA + col)*nRowA + (startRowA + row)];
    }
  }

}

// Convert a matrix to lower triangular
void mkLT(double *A, int n){
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < i; ++j){
      A[i * n + j] = 0.0;
    }
  }
}

// Solve linear system with upper-triangular Cholesky
void mysolveUT(double *A, double *b, int n){

  int info = 0;
  char const *upper = "U";
  char const *trans = "T";
  char const *ntrans = "N";
  char const *nunit = "N";
  int incx = 1;     // Increment for x

  F77_NAME(dpotrf)(upper, &n, A, &n, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}
  F77_NAME(dtrsv)(upper, trans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);
  F77_NAME(dtrsv)(upper, ntrans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);

}

// Solve linear system with lower-triangular Cholesky
void mysolveLT(double *A, double *b, int n){

  int info = 0;
  char const *lower = "L";
  char const *trans = "T";
  char const *ntrans = "N";
  char const *nunit = "N";
  int incx = 1;     // Increment for x

  F77_NAME(dpotrf)(lower, &n, A, &n, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}
  F77_NAME(dtrsv)(lower, ntrans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, trans, nunit, &n, A, &n, b, &incx FCONE FCONE FCONE);

}

// Print a matrix with entry type double
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

// Print a vector with entry type double
void printVec(double *m, int n){

  Rprintf("\t");
  for(int j = 0; j < n; j++){
    Rprintf("%.3f\t", m[j]);
  }
  Rprintf("\n");
}

// Print a vector with entry type integer
void printVec(int *m, int n){

  Rprintf(" ");
  for(int j = 0; j < n; j++){
    Rprintf("%i ", m[j]);
  }
  Rprintf("\n");
}

// Create lower-triangular spatial correlation matrix
void spCorLT(double *D, int n, double *theta, std::string &corfn, double *C){
  int i,j;

  if(corfn == "exponential"){

    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){
        C[i*n + j] = theta[0] * exp(-1.0 * theta[1] * D[i*n + j]);
      }
    }

  }else if(corfn == "matern"){

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
    perror("c++ error: corfn is not correctly specified");
  }
}

// Create full spatial correlation matrix
void spCorFull(double *D, int n, double *theta, std::string &corfn, double *C){
  int i,j;

  if(corfn == "exponential"){

    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){
        C[i*n + j] = theta[0] * exp(-1.0 * theta[1] * D[i*n + j]);
        C[j*n + i] = C[i*n + j];
      }
    }

  }else if(corfn == "matern"){

    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){
        if(D[i*n + j] * theta[0] > 0.0){
          C[i*n + j] = pow(D[i*n + j] * theta[0], theta[1]) / (pow(2, theta[1] - 1) * gammafn(theta[1])) * bessel_k(D[i*n + j] * theta[0], theta[1], 1.0);
          C[j*n + i] = C[i*n + j];
        }else{
          C[i*n + j] = 1.0;
        }
      }
    }

  }else{
    perror("c++ error: corfn is not correctly specified");
  }
}

// Fill a double vector with zeros
void zeros(double *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0.0;
}

// Fill an integer vector with zeros
void zeros(int *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0;
}
