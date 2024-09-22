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

// Copy a matrix excluding the i-th row
void copyMatrixDelRow(double *M1, int nRowM1, int nColM1, double *M2, int exclude_index){

  int i = 0, j = 0, new_index = 0;

  if(exclude_index < 0 || exclude_index > nRowM1){
    perror("Row index to exclude is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      for(i = 0; i < nRowM1; i++){
        if(i == exclude_index) continue;
        M2[new_index++] = M1[j*nRowM1 + i];
      }
    }
  }
}

// Copy a matrix excluding a row block
void copyMatrixDelRowBlock(double *M1, int nRowM1, int nColM1, double *M2, int exclude_start, int exclude_end){

  int i = 0, j = 0, new_index = 0;

  if(exclude_start > exclude_end || exclude_start == exclude_end){
    perror("Start index must be at least 1 less than End index.");
  }

  if(exclude_start < 0 || exclude_end > nRowM1){
    perror("Row index to exclude is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      for(i = 0; i < nRowM1; i++){
        if(i < exclude_start || i > exclude_end){
          M2[new_index++] = M1[j*nRowM1 + i];
        }
      }
    }
  }
}

// Copy a matrix deleting ith row and jth column
void copyMatrixDelRowCol(double *M1, int nRowM1, int nColM1, double *M2, int del_indexRow, int del_indexCol){

  int i = 0, j = 0, new_index = 0;

  if(del_indexRow < 0 || del_indexRow > nRowM1){
    perror("Row index to delete is out of bounds.");
  }else if(del_indexCol < 0 || del_indexCol > nColM1){
    perror("Column index to delete is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      if(j == del_indexCol) continue;
      for(i = 0; i < nRowM1; i++){
        if(i == del_indexRow) continue;
        M2[new_index++] = M1[j*nRowM1 + i];
      }
    }
  }
}

// Copy a matrix deleting ith row and jth column
void copyMatrixDelRowColBlock(double *M1, int nRowM1, int nColM1, double *M2,
                              int delRow_start, int delRow_end, int delCol_start, int delCol_end){

  int i = 0, j = 0, new_index = 0;

  if(delRow_start > delRow_end || delRow_start == delRow_end){
    perror("Row Start index must be at least 1 less than End index.");
  }

    if(delCol_start > delCol_end || delCol_start == delCol_end){
    perror("Column Start index must be at least 1 less than End index.");
  }

  if(delRow_start < 0 || delRow_end > nRowM1){
    perror("Row indices to delete are out of bounds.");
  }else if(delCol_start < 0 || delCol_end > nColM1){
    perror("Column indices to delete is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      if(j < delCol_start || j > delCol_end){
        for(i = 0; i < nRowM1; i++){
          if(i < delRow_start || i > delRow_end){
            M2[new_index++] = M1[j*nRowM1 + i];
          }
        }
      }
    }
  }
}

// Copy a block of rows of a matrix to another matrix
void copyMatrixRowBlock(double *M1, int nRowM1, int nColM1, double *M2, int copy_start, int copy_end){

  int i = 0, j = 0, new_index = 0;

  if(copy_start > copy_end || copy_start == copy_end){
    perror("Start index must be at least 1 less than End index.");
  }

  if(copy_start < 0 || copy_end > nRowM1){
    perror("Row indices to copy is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      for(i = 0; i < nRowM1; i++){
        if(i > copy_start - 1 && i < copy_end + 1){
          M2[new_index++] = M1[j*nRowM1 + i];
        }
      }
    }
  }

}

// Copy a row of a matrix to a vector
void copyMatrixRowToVec(double *M, int nRowM, int nColM, double *vec, int copy_index){

  int j = 0;

  // if(copy_index < 0 || copy_index > nRowM){
  //   perror("Row index to copy is out of bounds.");
  // }else{

  // }

  for(j = 0; j < nColM; j++){
    vec[j] = M[nRowM*j + copy_index];
  }

}

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

// Copy a vector excluding a block with start and end indices
void copyVecBlock(double *v1, double *v2, int n, int copy_start, int copy_end){

  int i = 0, j = 0;

  if(copy_start > copy_end || copy_start == copy_end){
    perror("Start index must be at least 1 less than End index.");
  }
  if(copy_start < 0 || copy_end > n){
    perror("Index to delete is out of bounds.");
  }else{
    for(i = 0; i < n; i++){
      if(i > copy_start - 1 && i < copy_end + 1){
        v2[j++] = v1[i];
      }
    }
  }
}

// Copy a vector excluding a block with start and end indices
void copyVecExcludingBlock(double *v1, double *v2, int n, int exclude_start, int exclude_end){

  int i = 0, j = 0;

  if(exclude_start > exclude_end || exclude_start == exclude_end){
    perror("Start index must be at least 1 less than End index.");
  }
  if(exclude_start < 0 || exclude_end > n){
    perror("Index to delete is out of bounds.");
  }else{
    for(i = 0; i < n; i++){
      if(i < exclude_start || i > exclude_end){
        v2[j++] = v1[i];
      }
    }
  }
}

// Copy a vector excluding the i-th entry
void copyVecExcludingOne(double *v1, double *v2, int n, int exclude_index){

  int i = 0, j = 0;

  if(exclude_index < 0 || exclude_index > n){
    perror("Index to delete is out of bounds.");
  }else{
    for(i = 0; i < n; i++){
      if(i != exclude_index){
        v2[j++] = v1[i];
      }
    }
  }
}

// Find maximum element in an integer vector
int findMax(int *a, int n){

  int i;
  int a_max = a[0];
  for(i = 1; i < n; i++){
    if(a[i] > a_max){
      a_max = a[i];
    }
  }

  return a_max;
}

// Find minimum element in an integer vector
int findMin(int *a, int n){

  int i;
  int a_min = a[0];
  for(i = 1; i < n; i++){
    if(a[i] < a_min){
      a_min = a[i];
    }
  }

  return a_min;
}

// Function to compute inverse-logit function
double inverse_logit(double x){
  return 1.0 / (1.0 + exp(-x));
}

// Function to compute log(x/(1-x)) for a given x
double logit(double x){
  return log(x) - log(1.0 - x);
}

// Function to compute logMeanExp of a vector
double logMeanExp(double *a, int n){

  int i;

  if(n == 0){
    perror("Vector of log values have 0 length.");
  }

  // Find maximum value in input vector
  double a_max = a[0];
  for(i = 1; i < n; i++){
    if(a[i] > a_max){
      a_max = a[i];
    }
  }

  // Find sum of adjusted exponentials; sum(exp(a_i - a_max))
  double sum_adj = 0.0;
  for(i = 0; i < n; i++){
    sum_adj += exp(a[i] - a_max);
  }

  // Find log-mean-exp; log(sum(exp(a_i))) - log(n)
  return a_max + log(sum_adj) - log(n);

}

// make partition for K-fold cross-validation, return partition start and end indices
void mkCVpartition(int n, int K, int *start_vec, int *end_vec, int *size_vec){

  int base_size = 0;             // Base-size of each partition
  int remainder = 0;             // Remaining elements to distribute
  int i, start = 0, end = 0;

  base_size = n / K;
  remainder = n % K;

  for(i = 0; i < K; i++){

    end = start + base_size - 1;

    size_vec[i] = base_size;

    if(remainder > 0){
      end++;
      remainder--;
      size_vec[i]++;
    }

    start_vec[i] = start;
    end_vec[i] = end;

    start = end + 1;
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
    Rprintf("%.7f\t", m[j]);
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
        C[i*n + j] = theta[0] * exp(-1.0 * theta[0] * D[i*n + j]);
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
