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

  if(exclude_start > exclude_end){
    perror("Exclude Start index must not exceed End index.");
  }

  if(include_start > include_end){
    perror("Copy Start index must not exceed End index.");
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

  if(exclude_start > exclude_end){
    perror("Start index must not exceed End index.");
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

// Copy a matrix excluding a row block in every n-th row block
void copyMatrixDelRowBlock_vc(double *M1, int nRowM1, int nColM1, double *M2, int exclude_start, int exclude_end, int rep){

  int i = 0, j = 0, new_index = 0;

  if(exclude_start > exclude_end){
    perror("Start index must not exceed End index.");
  }

  if(exclude_start < 0 || exclude_end > nRowM1*rep){
    perror("Row index to exclude is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      for(i = 0; i < nRowM1; i++){
        if(i % rep < exclude_start || i % rep > exclude_end){
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

  if(delRow_start > delRow_end){
    perror("Row Start index must not exceed End index.");
  }

    if(delCol_start > delCol_end){
    perror("Column Start index must not exceed End index.");
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

// Copy a matrix deleting a row and column block within each repxrep block
void copyMatrixDelRowColBlock_vc(double *M1, int nRowM1, int nColM1, double *M2, int delRow_start, int delRow_end,
                                 int delCol_start, int delCol_end, int rep){

  int i = 0, j = 0, new_index = 0;

  if(delRow_start > delRow_end){
    perror("Row Start index must not exceed End index.");
  }

  if(delCol_start > delCol_end){
    perror("Column Start index must not exceed End index.");
  }

  if(delRow_start < 0 || delRow_end > nRowM1){
    perror("Row indices to delete are out of bounds.");
  }else if(delCol_start < 0 || delCol_end > nColM1){
    perror("Column indices to delete is out of bounds.");
  }else{
    for(j = 0; j < nColM1; j++){
      if(j % rep < delCol_start || j % rep > delCol_end){
        for(i = 0; i < nRowM1; i++){
          if(i % rep < delRow_start || i % rep > delRow_end){
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

  if(copy_start > copy_end){
    perror("Start index must not exceed End index.");
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

  if(copyCol_start > copyCol_end){
    perror("Column Start index must not exceed End index.");
  }

  if(copyRow_start > copyRow_end){
    perror("Row Start index must not exceed End index.");
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

  if(copy_start > copy_end){
    perror("Start index must not exceed End index.");
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

  if(exclude_start > exclude_end){
    perror("Start index must not exceed End index.");
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

// Draw log(G), G ~ Gamma(shape, 1), without forming G when it could underflow. For shape >= 1 this is
// log(rgamma(shape, 1)) (G cannot underflow). For shape < 1 it uses the identity: if G1 ~ Gamma(shape + 1, 1)
// and U ~ Uniform(0, 1) independently, then G1*U^(1/shape) ~ Gamma(shape, 1); on the log scale
//   log(G) = log(G1) + log(U)/shape,
// which is finite for any shape > 0 (unif_rand() lies in the open interval (0, 1)). A direct
// rgamma(shape, 1) underflows to 0 with probability about 2^(-1074*shape)/Gamma(shape + 1), e.g. ~6e-4 at
// shape = 0.01. Exact in distribution; uses one extra uniform draw when shape < 1.
double rlogGamma(double shape){
  if(shape >= 1.0){
    return log(rgamma(shape, 1.0));
  }
  double lg1 = log(rgamma(shape + 1.0, 1.0));
  return lg1 + log(unif_rand()) / shape;
}

// Draw logit(B), B ~ Beta(a, b), as log(G1) - log(G2) with G1 ~ Gamma(a, 1) and G2 ~ Gamma(b, 1)
// independent (B = G1/(G1 + G2), so B/(1 - B) = G1/G2). Unlike logit(rbeta(a, b)), this never forms
// 1 - B, which rounds to 0 (B rounds to 1) with probability about pbeta(2^-53, b, a), e.g. ~5e-7 at
// (a, b) = (1.4, 0.4) and ~2e-5 at (1.3, 0.3); nor B itself, which can underflow when a is small.
// Exact in distribution.
double rlogitBeta(double a, double b){
  double lg1 = rlogGamma(a);
  double lg2 = rlogGamma(b);
  return lg1 - lg2;
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

// Print a matrix with entry type double
void printMtrx(double *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("% .2f\t", m[j*nRow+i]);
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

// Spatial correlation kernels. The kernel type and its constants are resolved once per
// matrix, not once per pair:
//   exponential: exp(-phi*d)
//   matern:      (phi*d)^nu / (2^(nu-1)*Gamma(nu)) * K_nu(phi*d), with the exact closed forms
//                nu = 0.5: exp(-x); nu = 1.5: (1 + x)*exp(-x); nu = 2.5: (1 + x + x^2/3)*exp(-x),
//                x = phi*d; otherwise K_nu is evaluated by bessel_k_ex with a preallocated work array.
// kernel codes
#define SPCOR_EXPONENTIAL 0
#define SPCOR_MATERN_05 1
#define SPCOR_MATERN_15 2
#define SPCOR_MATERN_25 3
#define SPCOR_MATERN 4

static int spCorCode(std::string &corfn, double nu){
  if(corfn == "exponential"){
    return SPCOR_EXPONENTIAL;
  }else if(corfn == "matern"){
    if(nu == 0.5) return SPCOR_MATERN_05;
    if(nu == 1.5) return SPCOR_MATERN_15;
    if(nu == 2.5) return SPCOR_MATERN_25;
    return SPCOR_MATERN;
  }else{
    Rf_error("c++ error: corfn is not correctly specified");
  }
  return -1;
}

// correlation at distance d; cnst = 1/(2^(nu-1)*Gamma(nu)), bk = work array of length floor(nu)+1
static inline double spCorEval(double d, int code, double phi, double nu, double cnst, double *bk){
  double x = phi * d;
  switch(code){
  case SPCOR_EXPONENTIAL:
  case SPCOR_MATERN_05:
    return exp(-x);
  case SPCOR_MATERN_15:
    return (1.0 + x) * exp(-x);
  case SPCOR_MATERN_25:
    return (1.0 + x + x * x / 3.0) * exp(-x);
  default:
    if(x > 0.0){
      return cnst * pow(x, nu) * bessel_k_ex(x, nu, 1.0, bk);
    }else{
      return 1.0;
    }
  }
}

// constants and work array of the general Matern kernel (no-op for the other kernels)
static double *spCorSetup(int code, double nu, double *cnst){
  double *bk = NULL;
  *cnst = 0.0;
  if(code == SPCOR_MATERN){
    *cnst = 1.0 / (pow(2.0, nu - 1.0) * gammafn(nu));
    bk = (double *) R_alloc((size_t) floor(nu) + 1, sizeof(double));
  }
  return bk;
}

// Euclidean distance between row i of the n x p matrix A and row j of the m x p matrix B
static inline double spDist(double *A, int n, int i, double *B, int m, int j, int p){
  int k;
  double dist = 0.0, dtemp = 0.0;
  for(k = 0; k < p; k++){
    dtemp = A[k * n + i] - B[k * m + j];
    dist += dtemp * dtemp;
  }
  return sqrt(dist);
}

// Create nxn full spatial correlation matrix from n x p coordinates
void spCorFull2(int n, int p, double *coords_sp, double *theta, std::string &corfn, double *C){
  int i, j;
  double nu = (corfn == "matern") ? theta[1] : 0.0;
  int code = spCorCode(corfn, nu);
  double cnst = 0.0;
  double *bk = spCorSetup(code, nu, &cnst);

  for(i = 0; i < n; i++){
    C[i*n + i] = 1.0;
    for(j = i + 1; j < n; j++){
      C[i*n + j] = spCorEval(spDist(coords_sp, n, i, coords_sp, n, j, p), code, theta[0], nu, cnst, bk);
      C[j*n + i] = C[i*n + j];
    }
  }
}

// Create full spatial-temporal correlation matrix
void sptCorFull(int n, int p, double *coords_sp, double *coords_tm, double *theta, std::string &corfn, double *C){
  int i, j, k;
  double sp_dist, tm_dist;

  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      sp_dist = 0.0;
      tm_dist = 0.0;

      // find spatial distance
      for(k = 0; k < p; k++){
        sp_dist += pow(coords_sp[k * n + i] - coords_sp[k * n + j], 2);
      }
      sp_dist = sqrt(sp_dist);

      // find temporal distance
      tm_dist = pow(coords_tm[i] - coords_tm[j], 2);
      tm_dist = sqrt(tm_dist);

      // evaluate correlation kernel
      if(corfn == "gneiting-decay"){
        C[i * n + j] = gneiting_spt_decay(sp_dist, tm_dist, theta[0], theta[1]);
        C[j * n + i] = C[i * n + j];
      }else{
        perror("c++ error: corfn is incorrectly specified");
      }

    }
  }
}

// Create nxn' spatial cross-correlation matrix
void spCorCross(int n, int n_prime, int p, double *coords_sp, double *coords_sp_prime, double *theta, std::string &corfn, double *C){
  int i, j;
  double nu = (corfn == "matern") ? theta[1] : 0.0;
  int code = spCorCode(corfn, nu);
  double cnst = 0.0;
  double *bk = spCorSetup(code, nu, &cnst);

  for(j = 0; j < n_prime; j++){
    for(i = 0; i < n; i++){
      C[j * n + i] = spCorEval(spDist(coords_sp, n, i, coords_sp_prime, n_prime, j, p), code, theta[0], nu, cnst, bk);
    }
  }
}

// Create nxn' cross-correlation spatial-temporal matrix
void sptCorCross(int n, int n_prime, int p, double *coords_sp, double *coords_tm, double *coords_sp_prime, double *coords_tm_prime, double *theta, std::string &corfn, double *C){
  int i, j, k;
  double sp_dist, tm_dist;

  for(i = 0; i < n; i++){
    for(j = 0; j < n_prime; j++){
      sp_dist = 0.0;
      tm_dist = 0.0;

      // find spatial distance
      for(k = 0; k < p; k++){
        sp_dist += pow(coords_sp[k * n + i] - coords_sp_prime[k * n_prime + j], 2);
      }
      sp_dist = sqrt(sp_dist);

      // find temporal distance
      tm_dist = pow(coords_tm[i] - coords_tm_prime[j], 2);
      tm_dist = sqrt(tm_dist);

      // evaluate correlation kernel
      if(corfn == "gneiting-decay"){
        C[j * n + i] = gneiting_spt_decay(sp_dist, tm_dist, theta[0], theta[1]);
      }else{
        perror("c++ error: corfn is incorrectly specified");
      }

    }
  }
}


// gneiting-decay spatio-temporal correlation function (Gneiting and Guttorp 2010)
double gneiting_spt_decay(double dist_s, double dist_t, double phi_s, double phi_t){

  double dist_t_sq = pow(dist_t, 2);
  double tmp = 0.0;
  tmp = (phi_t * dist_t_sq) + 1.0;

  return (1.0 / tmp) * exp(- (phi_s * dist_s) / sqrt(tmp));

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

// WARNING: the following function has the transpose case erroneous
// Function for sparse matrix-vector multiplication for varying coefficients models
void lmulv_XTilde_VC(const char *trans, int n, int r, double *XTilde, double *v, double *res){

  int i = 0, j = 0;
  const int inc_n = n;

  if(strcmp(trans, "N") == 0){
    for(i = 0; i < n; i++){
      res[i] = F77_CALL(ddot)(&r, &XTilde[i], &inc_n, &v[i], &inc_n);
    }
  }else if(strcmp(trans, "T") == 0){
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        res[i*n + j] = XTilde[i*n + j] * v[j];
      }
    }
  }else{
    perror("lmulv_XTilde_VC: Invalid transpose argument.");
  }

}

// Function for sparse matrix-matrix multiplication for varying coefficients models
void lmulm_XTilde_VC(const char *trans, int n, int r, int k, double *XTilde, double *A, double *res){

  int i = 0, j = 0, l = 0;
  const int inc_n = n;

  if(strcmp(trans, "N") == 0){
    for(i = 0; i < n; i++){
      for(j = 0; j < k; j++){
        res[j * n + i] = F77_CALL(ddot)(&r, &XTilde[i], &inc_n, &A[j * n * r + i], &inc_n);
      }
    }
  }else if(strcmp(trans, "T") == 0){
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        for(l = 0; l < k; l++){
          res[n*r*l + i*n + j] = XTilde[i*n + j] * A[n*l + j];
        }
      }
    }
  }else{
    perror("lmulm_XTilde_VC: Invalid transpose argument.");
  }

}

void rmul_Vz_XTildeT(int n, int r, double *XTilde, double *Vz, double *res, std::string &processtype){

  int i = 0, j = 0, l = 0;

  if(processtype == "independent.shared" || processtype == "multivariate"){
    for(l = 0; l < r; l++){
      for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
          res[n*r*j + l*n + i] = Vz[j*n + i] * XTilde[l*n + j];
        }
      }
    }
  }else if(processtype == "independent"){
    for(l = 0; l < r; l++){
      for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
          res[n*r*j + l*n + i] = Vz[n*n*l + j*n + i] * XTilde[l*n + j];
        }
      }
    }
  }

}

// Function for drawing a sample from a Inverse-Wishart distribution
// Bartlett factor of a Wishart W_r(nu, I) draw (Bartlett 1939; rectangular coordinates of Mahalanobis, Bose and
// Roy 1937): A lower triangular with A[i,i] = sqrt(chi-square(nu - i)) and standard normal entries below the
// diagonal, so that A*t(A) ~ W_r(nu, I). The random numbers consumed depend on r and nu only.
void rWishartBartlett(int r, double nu, double *A){

  int i = 0, j = 0;
  int rr = r * r;

  zeros(A, rr);
  // Fill diagonal with chi-square distributed values
  for(i = 0; i < r; i++){
    A[i * r + i] = sqrt(rchisq(nu-i));   //sqrt(rgamma(0.5*(nu - i), 2.0));
  }

  // Fill lower triangle with standard normal variates
  for(i = 1; i < r; i++){
    for(j = 0; j < i; j++){
      A[j * r + i] = rnorm(0.0, 1.0);
    }
  }

}

// Inverse-Wishart draw from a Bartlett factor A (see rWishartBartlett):
// Sigma = inv(L*A*t(A)*t(L)), L = cholinvIWscale (lower; its upper triangle is zeroed). Deterministic; tmp_rr is
// r x r workspace and may be the same array as A (A is not used after its first product).
// Returns 0 on success, 1 if the factorization or inversion failed (Sigma is then incomplete).
int invWishartFromBartlett(int r, double *A, double *cholinvIWscale, double *Sigma, double *tmp_rr){

  int info = 0;
  int i = 0, j = 0;
  char const *lower = "L";
  char const *ntran = "N";
  char const *ytran = "T";
  char const *rside = "R";
  const double one = 1.0;
  const double zero = 0.0;

  // Sigma = A * t(A)
  F77_NAME(dsyrk)(lower, ntran, &r, &r, &one, A, &r, &zero, Sigma, &r FCONE FCONE);
  // Sigma = cholinvIWscale * Sigma * t(cholinvIWscale)
  mkLT(cholinvIWscale, r);
  F77_NAME(dsymm)(rside, lower, &r, &r, &one, Sigma, &r, cholinvIWscale, &r, &zero, tmp_rr, &r FCONE FCONE);
  F77_NAME(dgemm)(ntran, ytran, &r, &r, &r, &one, tmp_rr, &r, cholinvIWscale, &r, &zero, Sigma, &r FCONE FCONE);

  F77_NAME(dpotrf)(lower, &r, Sigma, &r, &info FCONE); if(info != 0){return 1;}
  F77_NAME(dpotri)(lower, &r, Sigma, &r, &info FCONE); if(info != 0){return 1;}

  // make Sigma symmetric
  for(i = 1; i < r; i++){
    for(j = 0; j < i; j++){
      Sigma[i * r + j] = Sigma[j * r + i];
    }
  }

  return 0;

}

// Draw from an inverse-Wishart distribution (Bartlett factor, then the deterministic transform).
// Returns 0 on success, 1 on failure (see invWishartFromBartlett).
int rInvWishart(int r, double nu, double *cholinvIWscale, double *Sigma, double *tmp_rr){

  rWishartBartlett(r, nu, tmp_rr);
  return invWishartFromBartlett(r, tmp_rr, cholinvIWscale, Sigma, tmp_rr);

}

// Fit diagnostics (stored in the fit, reported by the R wrappers when verbose = TRUE; no extra factorization).
//
// Smallest and largest off-diagonal entry of the correlation matrix held in the lower triangle of the n x n
// column-major array A: the correlations of the two farthest-apart and of the two closest locations. One pass over
// the n(n-1)/2 entries below the diagonal; NA if n < 2.
void corOffDiagRange(double *A, int n, double *minCor, double *maxCor){

  int i, j;
  double lo = R_PosInf, hi = R_NegInf, a = 0.0;

  for(j = 0; j < n - 1; j++){
    for(i = j + 1; i < n; i++){
      a = A[(size_t) j * n + i];
      if(a < lo){ lo = a; }
      if(a > hi){ hi = a; }
    }
  }
  *minCor = (n < 2) ? NA_REAL : lo;
  *maxCor = (n < 2) ? NA_REAL : hi;

}

// Smallest relative Cholesky pivot min_i L[i,i]^2 / d[i] of the lower factor L (n x n, column-major) of a matrix
// with diagonal d (d[i] = dconst for all i if d is NULL). For a correlation-type matrix, L[i,i]^2 / d[i] is the part
// of the variance of variable i not explained by variables 1, ..., i-1; its inverse bounds the condition number below.
double minRelPivot(double *L, int n, double *d, double dconst){

  int i;
  double m = R_PosInf, piv = 0.0;

  for(i = 0; i < n; i++){
    piv = L[(size_t) i * n + i];
    piv = piv * piv / ((d == NULL) ? dconst : d[i]);
    if(piv < m){ m = piv; }
  }

  return m;

}

// Returns a copy of the named list list_r with the element "diagnostics" = c(min.pivot, min.cor, max.cor) appended.
SEXP appendDiagnostics(SEXP list_r, double minPivot, double minCor, double maxCor){

  int k, len = Rf_length(list_r);
  SEXP names_r = Rf_getAttrib(list_r, R_NamesSymbol);
  SEXP out_r = PROTECT(Rf_allocVector(VECSXP, len + 1));
  SEXP outNames_r = PROTECT(Rf_allocVector(STRSXP, len + 1));
  SEXP diag_r = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP diagNames_r = PROTECT(Rf_allocVector(STRSXP, 3));

  for(k = 0; k < len; k++){
    SET_VECTOR_ELT(out_r, k, VECTOR_ELT(list_r, k));
    SET_STRING_ELT(outNames_r, k, STRING_ELT(names_r, k));
  }
  REAL(diag_r)[0] = minPivot;
  REAL(diag_r)[1] = minCor;
  REAL(diag_r)[2] = maxCor;
  SET_STRING_ELT(diagNames_r, 0, Rf_mkChar("min.pivot"));
  SET_STRING_ELT(diagNames_r, 1, Rf_mkChar("min.cor"));
  SET_STRING_ELT(diagNames_r, 2, Rf_mkChar("max.cor"));
  Rf_setAttrib(diag_r, R_NamesSymbol, diagNames_r);
  SET_VECTOR_ELT(out_r, len, diag_r);
  SET_STRING_ELT(outNames_r, len, Rf_mkChar("diagnostics"));
  Rf_setAttrib(out_r, R_NamesSymbol, outNames_r);

  UNPROTECT(4);

  return out_r;

}

// As appendDiagnostics, for m correlation matrices (m processes): with m > 1, "diagnostics" is an m x 3 matrix with
// one row per process and columns min.pivot, min.cor, max.cor.
SEXP appendDiagnosticsRows(SEXP list_r, double *minPivot, double *minCor, double *maxCor, int m){

  if(m == 1){
    return appendDiagnostics(list_r, minPivot[0], minCor[0], maxCor[0]);
  }

  int k, len = Rf_length(list_r);
  SEXP names_r = Rf_getAttrib(list_r, R_NamesSymbol);
  SEXP out_r = PROTECT(Rf_allocVector(VECSXP, len + 1));
  SEXP outNames_r = PROTECT(Rf_allocVector(STRSXP, len + 1));
  SEXP diag_r = PROTECT(Rf_allocMatrix(REALSXP, m, 3));
  SEXP dimnames_r = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP colNames_r = PROTECT(Rf_allocVector(STRSXP, 3));

  for(k = 0; k < len; k++){
    SET_VECTOR_ELT(out_r, k, VECTOR_ELT(list_r, k));
    SET_STRING_ELT(outNames_r, k, STRING_ELT(names_r, k));
  }
  for(k = 0; k < m; k++){
    REAL(diag_r)[k] = minPivot[k];
    REAL(diag_r)[m + k] = minCor[k];
    REAL(diag_r)[2*m + k] = maxCor[k];
  }
  SET_STRING_ELT(colNames_r, 0, Rf_mkChar("min.pivot"));
  SET_STRING_ELT(colNames_r, 1, Rf_mkChar("min.cor"));
  SET_STRING_ELT(colNames_r, 2, Rf_mkChar("max.cor"));
  SET_VECTOR_ELT(dimnames_r, 1, colNames_r);
  Rf_setAttrib(diag_r, R_DimNamesSymbol, dimnames_r);
  SET_VECTOR_ELT(out_r, len, diag_r);
  SET_STRING_ELT(outNames_r, len, Rf_mkChar("diagnostics"));
  Rf_setAttrib(out_r, R_NamesSymbol, outNames_r);

  UNPROTECT(5);

  return out_r;

}

// User interrupt check that does not jump: R_CheckUserInterrupt() longjmps out of the C code, which would leak the
// R_chk_calloc buffers of the long leave-one-out / cross-validation loops. Run it inside R_ToplevelExec instead;
// returns 1 if the user asked to interrupt (the caller frees its memory and then stops with an error).
static void checkInterruptFn(void *dummy){
  R_CheckUserInterrupt();
}

int pendingInterrupt(){
  return !(R_ToplevelExec(checkInterruptFn, NULL));
}
