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

// Copy a matrix columns (submatrix) excluding a row block
void copyMatrixColDelRowBlock(double *M1, int nRowM1, int nColM1, double *M2,
                              int include_start, int include_end, int exclude_start, int exclude_end){

  int i = 0, j = 0, new_index = 0;

  if(exclude_start > exclude_end || exclude_start == exclude_end){
    perror("Exclude Start index must be at least 1 less than End index.");
  }

  if(include_start > exclude_start || include_start == include_end){
    perror("Copy Start index must be at least 1 less than End index.");
  }

  if(include_start < 0 || include_end > nColM1){
    perror("Column index to include is out of bounds.");
  }

  if(exclude_start < 0 || exclude_end > nRowM1){
    perror("Row index to exclude is out of bounds.");
  }else{
    for(j = include_start; j < include_end + 1; j++){
      for(i = 0; i < nRowM1; i++){
        if(i < exclude_start || i > exclude_end){
          M2[new_index++] = M1[j*nRowM1 + i];
        }
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


// Copy a matrix deleting a row and column block
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

// Copy a block (rows and columns) of a matrix to another matrix
void copyMatrixRowColBlock(double *M1, int nRowM1, int nColM1, double *M2,
                           int copyCol_start, int copyCol_end, int copyRow_start, int copyRow_end){

  int i = 0, j = 0, new_index = 0;

  if(copyCol_start > copyCol_end || copyCol_start == copyCol_end){
    perror("Column Start index must be at least 1 less than End index.");
  }

  if(copyRow_start > copyRow_end || copyRow_start == copyRow_end){
    perror("Row Start index must be at least 1 less than End index.");
  }

  if(copyRow_start < 0 || copyRow_end > nRowM1){
    perror("Row indices to copy is out of bounds.");
  }else if(copyCol_start < 0 || copyCol_end > nColM1){
    perror("Column indices to copy is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      if(j > copyCol_start - 1 && j < copyCol_end + 1){
        for(i = 0; i < nRowM1; i++){
          if(i > copyRow_start - 1 && i < copyRow_end + 1){
            M2[new_index++] = M1[j*nRowM1 + i];
          }
        }
      }
    }
  }

}

// Copy a column of a matrix to a vector
void copyMatrixColToVec(double *M, int nRowM, int nColM, double *vec, int copy_index){

  int i = 0;

  if(copy_index < 0 || copy_index > nColM){
    perror("Column index to copy is out of bounds.");
  }else{
    for(i = 0; i < nRowM; i++){
      vec[i] = M[nRowM*copy_index + i];
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

// Find maximum element in an double vector
double findMax(double *a, int n){

  int i;
  double a_max = a[0];
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

// Function to compute logSumExp of a vector
double logSumExp(double *a, int n){

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

  // Find log-sum-exp; log(sum(exp(a_i)))
  return a_max + log(sum_adj);

}

// Function to compute log-Weighted-Sum-Exp of a vector, weights in log-scale
double logWeightedSumExp(double *a, double *log_w, int n){

  int i;
  double log_num = 0, log_den = 0;

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
    sum_adj += exp(log_w[i] + (a[i] - a_max));
  }

  log_num = a_max + log(sum_adj);
  log_den = logSumExp(log_w, n);

  // Find log-sum-exp; log(sum(exp(a_i)))
  return log_num - log_den;

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
      Rprintf("%.2f\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }
}

// Print a vector with entry type double
void printVec(double *m, int n){

  Rprintf("\t");
  for(int j = 0; j < n; j++){
    Rprintf("%.2f\t", m[j]);
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

// Structure to hold both value and original index
typedef struct{
  double value;
  int index;
} IndexedValue;

// // Comparator function for qsort
// int compare(const void *a, const void *b){
//   double diff = ((IndexedValue *)a)->value - ((IndexedValue *)b)->value;
//   if(diff < 0) return -1;
//   if(diff > 0) return 1;
//   return 0;
// }

// Comparator function for qsort
int compare(const void *a, const void *b) {
    // Cast the pointers to IndexedValue and compare the values
    IndexedValue *ia = (IndexedValue *)a;
    IndexedValue *ib = (IndexedValue *)b;

    if (ia->value < ib->value) return -1;   // Return -1 if first value is smaller
    if (ia->value > ib->value) return 1;    // Return 1 if first value is larger
    return 0;                               // Return 0 if both are equal
}

// Pure C function to sort a vector and return the order (indices)
void sort_with_order(double *vec, int n, double *sorted_vec, int *order) {

  // Create an array of IndexedValue to hold both values and their original indices
  IndexedValue *arr = (IndexedValue *)malloc(n * sizeof(IndexedValue));
  if(arr == NULL){
    perror("Memory allocation failed");
  }

  // Populate the arr with the values and their original indices
  for(int i = 0; i < n; i++){
    arr[i].value = vec[i];
    arr[i].index = i;  // Store original index
  }

  // Sort the arr based on the value
  qsort(arr, n, sizeof(IndexedValue), compare);

  // After sorting, store the indices (order) in the output array
  for(int i = 0; i < n; i++){
    sorted_vec[i] = arr[i].value;
    order[i] = arr[i].index;  // Store the original indices
  }

  // Free the allocated memory
  free(arr);

}


// Fit generalized Pareto on raw importance ratios and return stabilized weights
void ParetoSmoothedIR(double *raw_IR, int M, int n_samples, double *sorted_IR, int *order_ind, double *stable_IR,
                      double *results, double *tailIR, double *exp_tail, double *stable_tail){

  int i = 0, ind = 0;
  double cutoff = 0.0, exp_cutoff = 0.0, tmp = 0.0;
  double max_raw_IR = 0.0;
  double k_hat = 0.0, sigma_hat = 0.0, p_tail = 0.0;

  max_raw_IR = findMax(raw_IR, n_samples);
  // Shift raw importance ratios for safe exponentiation
  for(i = 0; i < n_samples; i++){
    raw_IR[i] = raw_IR[i] - max_raw_IR;
  }

  // Sort raw importance ratios and store original indices
  zeros(order_ind, n_samples);
  sort_with_order(raw_IR, n_samples, sorted_IR, order_ind);

  for(i = 0; i < M; i++){
    tailIR[i] = sorted_IR[n_samples - M + i];
  }

  cutoff = sorted_IR[n_samples - M - 1];
  exp_cutoff = exp(cutoff);

  for(i = 0; i < M; i++){
    exp_tail[i] = exp(tailIR[i]) - exp_cutoff;
  }

  if(M > 5){

    fitGeneralParetoDist(exp_tail, M, 1, 30, results);
    k_hat = results[0];
    sigma_hat = results[1];

    for(i = 0; i < M; i++){
      p_tail = (i + 1 - 0.5) / M;
      tmp = qGPD(p_tail, k_hat, sigma_hat);
      exp_tail[i] = tmp + exp_cutoff;
      tmp = log(exp_tail[i]);
      stable_tail[i] = tmp;
    }

  }

  for(i = 0; i < M; i++){
    sorted_IR[n_samples - M + i] = tailIR[i];
  }

  for(i = 0; i < n_samples; i++){
    ind = order_ind[i];
    stable_IR[ind] = sorted_IR[i];
  }

  // truncate at max of raw wts (i.e., 0 since max has been subtracted)
  for(i = 0; i < n_samples; i++){
    if(stable_IR[i] > 0){
      stable_IR = 0;
    }
  }

  // shift back log-weights
  for(i = 0; i < n_samples; i++){
    stable_IR[i] += max_raw_IR;
  }

}

// Fit generalized Pareto distribution on a sorted sample
// Algorithm is based on Zhang, J., and Stephens, M. A. (2009). A new and efficient
// estimation method for the generalized Pareto distribution. Technometrics 51, 316-325.
// Closely follows gpdfit function of the R package "loo".
void fitGeneralParetoDist(double *x, int n, int wip, int min_grid_pts, double *result){

  int i = 0;
  int sqrt_n = 0, first_quart = 0;
  int m = 0;
  double xstar = 0.0, theta_hat = 0.0;
  double n_d = 0.0, m_d = 0.0;
  double prior = 3;
  const int incOne = 1;
  double k_hat = 0.0, sigma_hat = 0.0;

  n_d = n;
  sqrt_n = sqrt(n);
  m = min_grid_pts + sqrt_n;
  m_d = m;

  first_quart = 0.5 + (n_d / 4);
  xstar = x[first_quart - 1];

  double *theta = (double *) R_chk_calloc(m, sizeof(double)); zeros(theta, m);
  double *l_theta = (double *) R_chk_calloc(m, sizeof(double)); zeros(l_theta, m);
  double *w_theta = (double *) R_chk_calloc(m, sizeof(double)); zeros(w_theta, m);

  for(i = 0; i < m; i++){
    theta[i] = (1.0 / x[n-1]) + ((1 - sqrt(m_d / (i + 1 - 0.5))) / (prior * xstar));
    l_theta[i] = n * lx(theta[i], x, n);
  }

  for(i = 0; i < m; i++){
    w_theta[i] = exp(l_theta[i] - logSumExp(l_theta, m));
  }

  theta_hat = F77_CALL(ddot)(&m, theta, &incOne, w_theta, &incOne);

  for(i = 0; i < n; i++){
    k_hat += log1p(- theta_hat * x[i]);
  }
  k_hat = k_hat / n;
  sigma_hat = - k_hat / theta_hat;

  if(wip){
    k_hat = (k_hat * n) / (n + 10) + (0.5 * 10) / (n + 10);
  }

  result[0] = k_hat;
  // Rprintf("%.7f\n", k_hat);
  result[1] = sigma_hat;

  R_chk_free(theta);
  R_chk_free(l_theta);
  R_chk_free(w_theta);

}

// Function to evaluate profile log-likelihood of generalized Pareto
double lx(double b, double *x, int n){

  int i = 0;
  double sum = 0.0, k = 0.0;

  for(i = 0; i < n; i++){
    sum += log1p(- b * x[i]);
  }
  k = - sum / n;

  return log(b / k) + k - 1;

}

// Function to find quantiles of generalized pareto
double qGPD(double p, double k, double sigma){

  double out = 0.0;

  out = sigma * expm1(- k * log1p(- p)) / k;

  return out;

}
