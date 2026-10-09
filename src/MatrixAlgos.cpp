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

// Rank-1 update of Cholesky factor; chol(alpha*LLt + beta*vvt), as appearing in Krause and Igel (2015).
// REFERENCE:
// Oswin Krause and Christian Igel. 2015. A More Efficient Rank-one Covariance Matrix Update for Evolution Strategies.
// In Proceedings of the 2015 ACM Conference on Foundations of Genetic Algorithms XIII (FOGA '15). Association for
// Computing Machinery, New York, NY, USA, 129–136. https://doi.org/10.1145/2725494.2725496
void cholRankOneUpdate(int n, double *L1, double alpha, double beta, double *v, double *L2, double *w){

  int j, k;
  const int incOne = 1;
  const double sqrtalpha = sqrt(alpha);

  double b = 0.0, gamma = 0.0;
  double tmp1 = 0.0, tmp2 = 0.0;

  F77_NAME(dcopy)(&n, v, &incOne, w, &incOne);
  b = 1.0;

  for(j = 0; j < n; j++){

    tmp1 = pow(L1[j*n + j], 2);    // tmp1 = L[jj]^2
    tmp1 = alpha * tmp1;           // tmp1 = alpha*L[jj]^2
    tmp2 = pow(w[j], 2);           // tmp2 = w[j]^2
    tmp2 = beta * tmp2;            // tmp2 = beta*w[j]^2
    gamma = tmp1 * b;              // gamma = alpha*L[jj]^2*b
    gamma = gamma + tmp2;          // gamma = alpha*L[jj]^2*b + beta*w[j]^2
    tmp2 = tmp2 / b;               // tmp2 = (beta/b)*w[j]^2
    tmp1 = tmp1 + tmp2;            // tmp1 = alpha*L[jj]^2 + (beta/b)*w[j]^2
    tmp2 = sqrt(tmp1);             // tmp2 = sqrt(alpha*L[jj]^2 + (beta/b)*w[j]^2)

    L2[j*n +j] = tmp2;             // obtain L'[jj]

    if(j < n - 1){
      for(k = j + 1; k < n; k++){

        tmp1 = sqrtalpha * L1[j*n +k];   // tmp1 = sqrt(alpha)*L[kj]
        tmp1 = tmp1 / L1[j*n +j];        // tmp1 = sqrt(alpha)*L[kj]/L[jj]
        tmp2 = w[j] * tmp1;              // tmp2 = w[j]*(sqrt(alpha)*L[kj]/L[jj])
        w[k] = w[k] - tmp2;              // w[k] = w[k] - w[j]*(sqrt(alpha)*L[kj]/L[jj])

        tmp2 = beta * w[j];              // tmp2 = beta*w[j]
        tmp2 = tmp2 / gamma;             // tmp2 = (beta*w[j])/gamma
        tmp2 = tmp2 * w[k];              // tmp2 = (beta*w[j])*w[k]/gamma
        tmp2 = tmp2 + tmp1;              // tmp2 = sqrt(alpha)*L[kj]/L[jj] + (beta*w[j])*w[k]/gamma
        L2[j*n + k] = L2[j*n +j] * tmp2; // obtain L'[kj]

      }
    }

    tmp1 = pow(w[j], 2);           // tmp1 = w[j]^2
    tmp1 = beta * tmp1;            // tmp1 = beta*w[j]^2
    tmp2 = pow(L1[j*n +j], 2);     // tmp2 = L[jj]^2
    tmp2 = alpha * tmp2;           // tmp2 = alpha*L[jj]^2
    tmp1 = tmp1 / tmp2;            // tmp1 = beta*(w[j]^2/(alpha*L[jj]^2))
    b = b + tmp1;

  }

}

// Cholesky factor update after deletion of a row/column where rank-1 updates are
// carried out using Krause and Igel (2015).
void cholRowDelUpdate(int n, double *L, int del, double *L1, double *w){

  int j, k;
  const int n1 = n - 1;
  const int incOne = 1;

  int nk = 0;
  int indexLjj = 0, indexLkj= 0;
  double b = 0.0, gamma = 0.0;
  double tmp1 = 0.0, tmp2 = 0.0;

  if(del == n - 1){

    copySubmat(L, n, n, L1, n1, n1, 0, 0, 0, 0, n1, n1);
    mkLT(L1, n1);

  }else if(del == 0){

    nk = n - 1;
    int delPlusOne = del + 1;
    // w = (double *) R_chk_realloc(w, nk * sizeof(double));
    F77_NAME(dcopy)(&n1, &L[1], &incOne, w, &incOne);
    b = 1.0;

    for(j = 0; j < nk; j++){

      indexLjj = mapIndex(j, j, nk, nk, delPlusOne, delPlusOne, n);
      tmp1 = pow(L[indexLjj], 2);     // tmp1 = L[jj]^2
      gamma = tmp1 * b;               // gamma = L[jj]^2*b
      tmp2 = pow(w[j], 2);            // tmp2 = w[j]^2
      gamma = gamma + tmp2;           // gamma = L[jj]^2*b + w[j]^2
      tmp2 = tmp2 / b;                // tmp2 = w[j]^2/b
      tmp1 = tmp1 + tmp2;             // tmp1 = L[jj]^2 + w[j]^2/b
      tmp2 = sqrt(tmp1);              // tmp2 = sqrt(L[jj]^2 + w[j]^2/b)
      L1[j*nk + j] = tmp2;            // obtain L'[jj]

      if(j < nk - 1){
        for(k = j + 1; k < nk; k++){

          indexLkj = mapIndex(k, j, nk, nk, delPlusOne, delPlusOne, n);
          tmp1 = L[indexLkj] / L[indexLjj];   // tmp1 = L[kj]/L[jj]
          tmp2 = tmp1 * w[j];                 // tmp2 = w[j]*L[kj]/L[jj]
          w[k] = w[k] - tmp2;                 // w[k] = w[k] - w[j]*L[kj]/L[jj]

          tmp2 = w[j] * w[k];                 // tmp2 = w[j]*w[k]
          tmp2 = tmp2 / gamma;                // tmp2 = w[j]*w[k]/gamma
          tmp1 = tmp1 + tmp2;                 // tmp1 = L[kj]/L[jj] + w[j]*w[k]/gamma
          tmp2 = tmp1 * L1[j*nk + j];         // tmp1 = L'[jj]*L[kj]/L[jj] + L'[jj]*w[j]*w[k]/gamma
          L1[j*nk + k] = tmp2;                // obtain L'[kj]

        }

        tmp1 = pow(w[j], 2);          // tmp1 = w[j]^2
        tmp2 = pow(L[indexLjj], 2);   // tmp2 = L[jj]^2
        tmp1 = tmp1 / tmp2;           // tmp1 = w[j]^2/L[jj]^2
        b = b + tmp1;                 // b = b + w[j]^2/L[jj]^2

      }

      mkLT(L1, n1);

    }  // End rank-one update for first row/column deletion

  }else if(0 < del && del < n - 1){

    int delPlusOne = del + 1;
    int indexL1 = 0;

    nk = n - delPlusOne;

    copySubmat(L, n, n, L1, n1, n1, 0, 0, 0, 0, del, del);
    copySubmat(L, n, n, L1, n1, n1, delPlusOne, 0, del, 0, nk, del);

    // w = (double *) R_chk_realloc(w, nk * sizeof(double));
    F77_NAME(dcopy)(&nk, &L[del*n + delPlusOne], &incOne, w, &incOne);
    b = 1.0;

    for(j = 0; j < nk; j++){

      indexLjj = mapIndex(j, j, nk, nk, delPlusOne, delPlusOne, n);
      tmp1 = pow(L[indexLjj], 2);     // tmp1 = L[jj]^2
      gamma = tmp1 * b;               // gamma = L[jj]^2*b
      tmp2 = pow(w[j], 2);            // tmp2 = w[j]^2
      gamma = gamma + tmp2;           // gamma = L[jj]^2*b + w[j]^2
      tmp2 = tmp2 / b;                // tmp2 = w[j]^2/b
      tmp1 = tmp1 + tmp2;             // tmp1 = L[jj]^2 + w[j]^2/b
      tmp2 = sqrt(tmp1);              // tmp2 = sqrt(L[jj]^2 + w[j]^2/b)
      indexL1 = mapIndex(j, j, nk, nk, del, del, n1);
      L1[indexL1] = tmp2;            // obtain L'[jj]

      if(j < nk - 1){
        for(k = j + 1; k < nk; k++){

          indexLkj = mapIndex(k, j, nk, nk, delPlusOne, delPlusOne, n);
          tmp1 = L[indexLkj] / L[indexLjj];   // tmp1 = L[kj]/L[jj]
          tmp2 = tmp1 * w[j];                 // tmp2 = w[j]*L[kj]/L[jj]
          w[k] = w[k] - tmp2;                 // w[k] = w[k] - w[j]*L[kj]/L[jj]

          tmp2 = w[j] * w[k];                 // tmp2 = w[j]*w[k]
          tmp2 = tmp2 / gamma;                // tmp2 = w[j]*w[k]/gamma
          tmp1 = tmp1 + tmp2;                 // tmp1 = L[kj]/L[jj] + w[j]*w[k]/gamma
          indexL1 = mapIndex(j, j, nk, nk, del, del, n1);
          tmp2 = tmp1 * L1[indexL1];         // tmp1 = L'[jj]*L[kj]/L[jj] + L'[jj]*w[j]*w[k]/gamma
          indexL1 = mapIndex(k, j, nk, nk, del, del, n1);
          L1[indexL1] = tmp2;                 // obtain L'[kj]

        }

        tmp1 = pow(w[j], 2);          // tmp1 = w[j]^2
        tmp2 = pow(L[indexLjj], 2);   // tmp2 = L[jj]^2
        tmp1 = tmp1 / tmp2;           // tmp1 = w[j]^2/L[jj]^2
        b = b + tmp1;                 // b = b + w[j]^2/L[jj]^2

      }

    }

    mkLT(L1, n1);

  }else{
    perror("Row/column deletion index out of bounds.");
  }

}

// Cholesky factor update after deletion of a block
// using rank-1 updates as given in Krause and Igel (2015).
void cholBlockDelUpdate(int n, double *L, int del_start, int del_end, double *L1, double *tmpL1, double *w){

  int j, k;
  const int incOne = 1;
  int case_id = 0, nk = 0, nkk = 0, nMnk = 0, nMnknMnk = 0;
  int del = 0, delEndPlusOne = 0;
  double b = 0, gamma = 0, tmp1 = 0, tmp2 = 0;
  int indexLjj = 0, indexLkj= 0;

  // Error handling
  // indices are 0-based; a block of size 1 (del_start == del_end) is valid
  if(del_start > del_end){
    perror("Block Start index must not exceed End index.");
    return;
  }
  if(del_start < 0 || del_end > n - 1){
    perror("Block index to delete is out of bounds.");
    return;
  }

  // Step 1: Determine if deletion case is terminal or intermediate
  if(del_start > 0 && del_end == n - 1){
    case_id = 1;                           // Lowest block deletion
  }else if(del_start == 0 && del_end < n - 1){
    case_id = 2;                           // First block deletion
  }else{
    case_id = 3;
  }

  if(case_id == 1){

    nk = del_end - del_start + 1;
    nMnk = n - nk;
    copySubmat(L, n, n, L1, nMnk, nMnk, 0, 0, 0, 0, nMnk, nMnk);
    mkLT(L1, nMnk);

  }else if(case_id == 2){

    nk = del_end - del_start + 1;
    nMnk = n - nk;
    nMnknMnk = nMnk * nMnk;
    delEndPlusOne = del_end + 1;

    copySubmat(L, n, n, tmpL1, nMnk, nMnk, delEndPlusOne, delEndPlusOne, 0, 0, nMnk, nMnk);

    for(del = del_start; del < delEndPlusOne; del++){

      F77_NAME(dcopy)(&nMnk, &L[del*n + delEndPlusOne], &incOne, w, &incOne);

      b = 1.0;

      for(j = 0; j < nMnk; j++){

        tmp1 = pow(tmpL1[j*nMnk + j], 2);  // tmp1 = L[jj]^2
        gamma = tmp1 * b;                  // gamma = L[jj]^2*b
        tmp2 = pow(w[j], 2);               // tmp2 = w[j]^2
        gamma = gamma + tmp2;              // gamma = L[jj]^2*b + w[j]^2
        tmp2 = tmp2 / b;                   // tmp2 = w[j]^2/b
        tmp1 = tmp1 + tmp2;                // tmp1 = L[jj]^2 + w[j]^2/b
        tmp2 = sqrt(tmp1);                 // tmp2 = sqrt(L[jj]^2 + w[j]^2/b)
        L1[j*nMnk + j] = tmp2;             // obtain L'[jj]

        if(j < nMnk - 1){
          for(k = j + 1; k < nMnk; k++){

            tmp1 = tmpL1[j*nMnk + k] / tmpL1[j*nMnk + j];  // tmp1 = L[kj]/L[jj]
            tmp2 = tmp1 * w[j];                            // tmp2 = w[j]*L[kj]/L[jj]
            w[k] = w[k] - tmp2;                            // w[k] = w[k] - w[j]*L[kj]/L[jj]

            tmp2 = w[j] * w[k];                            // tmp2 = w[j]*w[k]
            tmp2 = tmp2 / gamma;                           // tmp2 = w[j]*w[k]/gamma
            tmp1 = tmp1 + tmp2;                            // tmp1 = L[kj]/L[jj] + w[j]*w[k]/gamma
            tmp2 = tmp1 * L1[j*nMnk + j];                  // tmp1 = L'[jj]*L[kj]/L[jj] + L'[jj]*w[j]*w[k]/gamma
            L1[j*nMnk + k] = tmp2;                         // obtain L'[kj]

          }
        }

        tmp1 = pow(w[j], 2);                 // tmp1 = w[j]^2
        tmp2 = pow(tmpL1[j*nMnk + j], 2);    // tmp2 = L[jj]^2
        tmp1 = tmp1 / tmp2;                  // tmp1 = w[j]^2/L[jj]^2
        b = b + tmp1;                        // b = b + w[j]^2/L[jj]^2

      }

      if(del < del_end){
        F77_NAME(dcopy)(&nMnknMnk, &L1[0], &incOne, tmpL1, &incOne);
      }

    }

    mkLT(L1, nMnk);

  }else if(case_id == 3){

    nk = del_end - del_start + 1;
    nMnk = n - nk;
    delEndPlusOne = del_end + 1;
    nkk = n - delEndPlusOne;

    copySubmat(L, n, n, tmpL1, nMnk, nMnk, delEndPlusOne, delEndPlusOne, del_start, del_start, nkk, nkk);

    for(del = del_start; del < delEndPlusOne; del++){

      F77_NAME(dcopy)(&nkk, &L[del*n + delEndPlusOne], &incOne, w, &incOne);

      b = 1.0;

      for(j = 0; j < nkk; j++){

        indexLjj = mapIndex(j, j, nkk, nkk, del_start, del_start, nMnk);
        tmp1 = pow(tmpL1[indexLjj], 2);    // tmp1 = L[jj]^2
        gamma = tmp1 * b;                  // gamma = L[jj]^2*b
        tmp2 = pow(w[j], 2);               // tmp2 = w[j]^2
        gamma = gamma + tmp2;              // gamma = L[jj]^2*b + w[j]^2
        tmp2 = tmp2 / b;                   // tmp2 = w[j]^2/b
        tmp1 = tmp1 + tmp2;                // tmp1 = L[jj]^2 + w[j]^2/b
        tmp2 = sqrt(tmp1);                 // tmp2 = sqrt(L[jj]^2 + w[j]^2/b)
        L1[indexLjj] = tmp2;               // obtain L'[jj]

        if(j < nkk - 1){
          for(k = j + 1; k < nkk; k++){

            indexLkj = mapIndex(k, j, nkk, nkk, del_start, del_start, nMnk);
            tmp1 = tmpL1[indexLkj] / tmpL1[indexLjj];      // tmp1 = L[kj]/L[jj]
            tmp2 = tmp1 * w[j];                            // tmp2 = w[j]*L[kj]/L[jj]
            w[k] = w[k] - tmp2;                            // w[k] = w[k] - w[j]*L[kj]/L[jj]

            tmp2 = w[j] * w[k];                            // tmp2 = w[j]*w[k]
            tmp2 = tmp2 / gamma;                           // tmp2 = w[j]*w[k]/gamma
            tmp1 = tmp1 + tmp2;                            // tmp1 = L[kj]/L[jj] + w[j]*w[k]/gamma
            tmp2 = tmp1 * L1[indexLjj];                    // tmp1 = L'[jj]*L[kj]/L[jj] + L'[jj]*w[j]*w[k]/gamma
            L1[indexLkj] = tmp2;                           // obtain L'[kj]

          }
        }

        tmp1 = pow(w[j], 2);                 // tmp1 = w[j]^2
        tmp2 = pow(tmpL1[indexLjj], 2);      // tmp2 = L[jj]^2
        tmp1 = tmp1 / tmp2;                  // tmp1 = w[j]^2/L[jj]^2
        b = b + tmp1;                        // b = b + w[j]^2/L[jj]^2

      }

      if(del < del_end){
        copySubmat(L1, nMnk, nMnk, tmpL1, nMnk, nMnk, del_start, del_start, del_start, del_start, nkk, nkk);
      }

    }

    copySubmat(L, n, n, L1, nMnk, nMnk, 0, 0, 0, 0, del_start, del_start);
    copySubmat(L, n, n, L1, nMnk, nMnk, delEndPlusOne, 0, del_start, 0, nkk, del_start);
    mkLT(L1, nMnk);

  }else{
    perror("cholBlockDelUpdate error: Invalid case.");
  }

}

// get the Schur complement of xi-cov submatrix for GLM case
// Pre-processing for projGLM(). With P = inv(Vz + I) (applied through cholVzPlusI = chol(Vz + I)):
//   inv(inv(Vz) + I) = I - P,                      so D1invB1 = inv(inv(Vz) + I)*X = X - W, with W = P*X,
//   Schur(A1) = t(X)*X + inv(Vbeta) - t(X)*inv(inv(Vz) + I)*X = t(X)*W + inv(Vbeta),
//   DinvB_np  = W*inv(Schur(A1))                   (n x p; the transpose of inv(Schur(A1))*t(W)),
//   Schur(A)  = I/sigmaSqxi + Q,  Q = P - W*inv(Schur(A1))*t(W) = inv(Vz + I + X*Vbeta*t(X)).
// Outputs: out_pp = chol(Schur(A1)), out_nn = chol(Schur(A)) (lower triangles), DinvB_np, D1invB1.
// tmp_np is n x p workspace.
void cholSchurGLM(double *X, int n, int p, double sigmaSqxi, double *VbetaInv, double *cholVzPlusI,
                  double *tmp_np, double *DinvB_np, double *out_pp, double *out_nn, double *D1invB1){

  int np = n * p;
  int pp = p * p;
  int nn = n * n;
  int i;

  int info = 0;
  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  char const *lside = "L";
  char const *rside = "R";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;
  const double sigmaSqxiInv = 1.0 / sigmaSqxi;

  // W = inv(Vz + I)*X
  F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);                                                              // tmp_np = X
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholVzPlusI, &n, tmp_np, &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholVzPlusI, &n, tmp_np, &n FCONE FCONE FCONE FCONE); // tmp_np = W = inv(Vz+I)*X

  // D1invB1 = inv(inv(Vz) + I)*X = X - W
  F77_NAME(dcopy)(&np, X, &incOne, D1invB1, &incOne);
  F77_NAME(daxpy)(&np, &negone, tmp_np, &incOne, D1invB1, &incOne);                                               // D1invB1 = X - W

  // chol(Schur(A1)), Schur(A1) = t(X)*W + inv(Vbeta)
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np, &n, &zero, out_pp, &p FCONE FCONE);             // out_pp = t(X)*W
  F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, out_pp, &incOne);                                                 // out_pp = t(X)*W + inv(Vbeta)
  F77_NAME(dpotrf)(lower, &p, out_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}      // out_pp = chol(Schur(A1))

  // DinvB_np = W*inv(Schur(A1)) = W*t(Linv)*Linv, L = chol(Schur(A1)); the intermediate T = W*t(Linv)
  // gives W*inv(Schur(A1))*t(W) = T*t(T)
  F77_NAME(dcopy)(&np, tmp_np, &incOne, DinvB_np, &incOne);                                                       // DinvB_np = W
  F77_NAME(dtrsm)(rside, lower, ytran, nunit, &n, &p, &one, out_pp, &p, DinvB_np, &n FCONE FCONE FCONE FCONE);    // DinvB_np = T = W*t(Linv)

  // Schur(A) = I/sigmaSqxi + inv(Vz + I) - T*t(T) (lower triangle)
  F77_NAME(dcopy)(&nn, cholVzPlusI, &incOne, out_nn, &incOne);
  F77_NAME(dpotri)(lower, &n, out_nn, &n, &info FCONE); if(info != 0){perror("c++ error: dpotri failed\n");}      // out_nn = inv(Vz + I)
  F77_NAME(dsyrk)(lower, ntran, &n, &p, &negone, DinvB_np, &n, &one, out_nn, &n FCONE FCONE);                     // out_nn = inv(Vz + I) - T*t(T) = Q
  for(i = 0; i < n; i++){
    out_nn[i * n + i] += sigmaSqxiInv;                                                                            // out_nn = I/sigmaSqxi + Q
  }

  F77_NAME(dtrsm)(rside, lower, ntran, nunit, &n, &p, &one, out_pp, &p, DinvB_np, &n FCONE FCONE FCONE FCONE);    // DinvB_np = W*inv(Schur(A1)) and RETURN

  // Find Cholesky factor of Schur complement
  F77_NAME(dpotrf)(lower, &n, out_nn, &n, &info FCONE); if(info != 0){perror("c++ error: Schur dpotrf failed\n");}

}

// No memory allocation inside 'hot' function
void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1, double *tmp_pp,
                 double *tmp_np1, double *out_p, double *out_n, int LOO){

  int pp = p * p;
  // int np = n * p;

  int info = 0;
  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  char const *lside = "L";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;

  const double deltasqInv = 1.0 / deltasq;
  const double negdeltasqInv = - 1.0 / deltasq;

  if(LOO){

    F77_NAME(dcopy)(&n, v2, &incOne, tmp_n1, &incOne);                                                         // tmp_n1 = v2 = J
    F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);
    F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);                   // tmp_n1 = VyInv*J
    F77_NAME(dscal)(&n, &deltasq, tmp_n1, &incOne);                                                            // tmp_n1 = deltasq*VyInv*J

  }else{

    F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, v2, &incOne, &zero, tmp_n1, &incOne FCONE);                   // tmp_n1 = Vz*v2
    F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);
    F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);                   // tmp_n1 = VyInv*Vz*v2
    F77_NAME(dscal)(&n, &deltasq, tmp_n1, &incOne);                                                            // tmp_n1 = deltasq*VyInv*Vz*v2

  }

  F77_NAME(dcopy)(&n, tmp_n1, &incOne, out_n, &incOne);                                                        // out_n = tmp_n1 = inv(D)*v2
  F77_NAME(dcopy)(&p, v1, &incOne, tmp_p1, &incOne);                                                           // tmp_p1 = v1
  F77_NAME(dgemv)(ytran, &n, &p, &negdeltasqInv, X, &n, tmp_n1, &incOne, &one, tmp_p1, &incOne FCONE);         // tmp_p1 = v1 - t(B)*inv(D)*v2

  F77_NAME(dcopy)(&pp, VbetaInv, &incOne, tmp_pp, &incOne);                                                    // tmp_pp = VbetaInv
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &deltasqInv, X, &n, X, &n, &one, tmp_pp, &p FCONE FCONE);          // tmp_pp = A = (1/deltasq)*XtX+VbetaInv

  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &one, Vz, &n, X, &n, &zero, tmp_np1, &n FCONE FCONE);              // tmp_np1 = Vz*X = deltasq*Vz*B
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholVy, &n, tmp_np1, &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholVy, &n, tmp_np1, &n FCONE FCONE FCONE FCONE);  // tmp_np1 = deltasq*VyInv*Vz*B = inv(D)*B

  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &negdeltasqInv, X, &n, tmp_np1, &n, &one, tmp_pp, &p FCONE FCONE); // tmp_pp = Schur(A) = A - t(B)*inv(D)*B
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}   // chol(Schur(A))
  F77_NAME(dtrsv)(lower, ntran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                     // tmp_p1 = inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dcopy)(&p, tmp_p1, &incOne, out_p, &incOne);                                                        // out_p = first p elements of Mv

  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, tmp_p1, &incOne, &zero, tmp_n1, &incOne FCONE);                  // tmp_n1 = deltasq*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, tmp_n1, &incOne, &zero, tmp_n2, &incOne FCONE);                 // tmp_n2 = Vz * tmp_n1
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);                     // tmp_n2 = inv(D)*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(daxpy)(&n, &negone, tmp_n2, &incOne, out_n, &incOne);                                               // out_n = inv(D)*v2 - inv(D)*B*inv(Schur(A))*(v1-BtDinvB)

}

// Map the index of the (i, j)-th entry of B to the corresponding index in A, where B is a submatrix of A.
int mapIndex(int i, int j, int nRowB, int nColB, int startRowB, int startColB, int nRowA){

  // Calculate the row and column indices of B[i,j] in A
  int rowA = startRowB + i;
  int colA = startColB + j;

  // Calculate the index in column-major order
  int indexA = rowA + colA * nRowA;

  return indexA;
}

// projection operator for GLM
// Projection step of the GCM sampler for the spatial GLM. With P = inv(Vz + I) (through cholVzPlusI):
//   inv(D1)*v = v - P*v, and the n x n block of inv(D)*B applied to v is inv(D1)*v - D1invB1*(t(DinvB_np)*v),
// so neither Vz nor that n x n block is needed here.
void projGLM(double *X, int n, int p, double *v_eta, double *v_xi, double *v_beta, double *v_z,
             double *cholpSchur, double *cholnSchur, double sigmaSqxi, double *Lbeta, double *Lz,
             double *cholVzPlusI, double *D1invB1, double *DinvBnp, double *tmp_n, double *tmp_p){

  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;
  const double sigmaxiInv = 1.0 / sqrt(sigmaSqxi);

  // Find components of t(H)*v, where (3n+p)x1 vector v = [v_eta, v_xi, v_beta, v_z]
  F77_NAME(dscal)(&n, &sigmaxiInv, v_xi, &incOne);                                                 // v_xi = v_xi/sigmasqxi
  F77_NAME(daxpy)(&n, &one, v_eta, &incOne, v_xi, &incOne);                                        // v_xi = v_eta + v_xi/sigmasqxi

  F77_NAME(dtrsv)(lower, ytran, nunit, &p, Lbeta, &p, v_beta, &incOne FCONE FCONE FCONE);          // v_beta = LbetatInv*v_beta
  F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, v_eta, &incOne, &one, v_beta, &incOne FCONE);        // v_beta = Xt*v_eta + LbetaInv*v_beta

  F77_NAME(dtrsv)(lower, ytran, nunit, &n, Lz, &n, v_z, &incOne FCONE FCONE FCONE);                // v_z = LztInv*v_z
  F77_NAME(daxpy)(&n, &one, v_eta, &incOne, v_z, &incOne);                                         // v_z = v_eta + LztInv*v_z

  // Find inv(D1)*v22 = v22 - inv(Vz+I)*v22
  F77_NAME(dcopy)(&n, v_z, &incOne, tmp_n, &incOne);                                               // tmp_n = v22
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);     // tmp_n = inv(Vz+I)*v22
  F77_NAME(dscal)(&n, &negone, tmp_n, &incOne);
  F77_NAME(daxpy)(&n, &one, v_z, &incOne, tmp_n, &incOne);                                         // tmp_n = D1inv*v22

  // Find (v21 - t(B1)*D1inv*v22)
  F77_NAME(dgemv)(ytran, &n, &p, &negone, D1invB1, &n, v_z, &incOne, &one, v_beta, &incOne FCONE); // v21 = v21 - B1t*D1Inv*v22

  // Find inv(schurA1)*(v21 - t(B1)*D1inv*v22)
  F77_NAME(dtrsv)(lower, ntran, nunit, &p, cholpSchur, &p, v_beta, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &p, cholpSchur, &p, v_beta, &incOne FCONE FCONE FCONE);     // v_beta = inv(sA1)*(v21 - t(B1)*D1inv*v22)

  F77_NAME(dcopy)(&p, v_beta, &incOne, tmp_p, &incOne);                                            // tmp_p = inv(sA1)*(v21 - t(B1)*D1inv*v22)
  F77_NAME(dgemv)(ntran, &n, &p, &one, D1invB1, &n, tmp_p, &incOne, &zero, v_z, &incOne FCONE);    // v_z = D1invB1*inv(sA1)*(v21 - t(B1)*D1inv*v22)
  F77_NAME(dscal)(&n, &negone, v_z, &incOne);                                                      // v_z = -D1invB1*inv(sA1)*(v21 - t(B1)*D1inv*v22)
  F77_NAME(daxpy)(&n, &one, tmp_n, &incOne, v_z, &incOne);                                         // v_z = D1inv*v22 - D1invB1*inv(sA1)*(v21 - t(B1)*D1inv*v22)

  // Find inv(schurA)*(v1 - BtDInvv2)
  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, v_beta, &incOne, &zero, tmp_n, &incOne FCONE);       // tmp_n = X*v_beta
  F77_NAME(daxpy)(&n, &one, v_z, &incOne, tmp_n, &incOne);                                         // tmp_n = X*v_beta + v_z
  F77_NAME(daxpy)(&n, &negone, tmp_n, &incOne, v_xi, &incOne);                                     // v_xi = v_xi - (X*v_beta + v_z)

  // Find v_xi
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholnSchur, &n, v_xi, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholnSchur, &n, v_xi, &incOne FCONE FCONE FCONE);       // v_xi = inv(sA)*(v1-BtDInvv2)

  // Find DInvB*inv(schurA1)*(v1 - BtDInvv2): p x 1 block t(DinvB_np)*v_xi; n x 1 block
  // inv(D1)*v_xi - D1invB1*(t(DinvB_np)*v_xi), with inv(D1)*v_xi = v_xi - inv(Vz+I)*v_xi
  F77_NAME(dgemv)(ytran, &n, &p, &one, DinvBnp, &n, v_xi, &incOne, &zero, tmp_p, &incOne FCONE);   // tmp_p = t(DinvB_np)*v_xi
  F77_NAME(dcopy)(&n, v_xi, &incOne, tmp_n, &incOne);
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);     // tmp_n = inv(Vz+I)*v_xi
  F77_NAME(dscal)(&n, &negone, tmp_n, &incOne);
  F77_NAME(daxpy)(&n, &one, v_xi, &incOne, tmp_n, &incOne);                                        // tmp_n = inv(D1)*v_xi
  F77_NAME(dgemv)(ntran, &n, &p, &negone, D1invB1, &n, tmp_p, &incOne, &one, tmp_n, &incOne FCONE); // tmp_n = inv(D1)*v_xi - D1invB1*tmp_p

  // Find v_beta, v_z
  F77_NAME(daxpy)(&p, &negone, tmp_p, &incOne, v_beta, &incOne);
  F77_NAME(daxpy)(&n, &negone, tmp_n, &incOne, v_z, &incOne);

}

// Function to transpose a matrix in column-major form

// Function to transpose a matrix in column-major form from upper-tri to lower-tri
void upperTri_lowerTri(double *M, int n){

  int i = 0, j = 0;

  for(j = 0; j < n; j++){
    for(i = 0; i < n; i++){
      if(i < j){
        M[i*n + j] = M[j*n + i];
      }
    }
  }
}

// Function for priming (pre-proprocessing) step for varying-coefficients model
// Pre-processing for projGLMvc(). With XTf = [diag(XTilde[,1]) ... diag(XTilde[,r])] (n x nr), Vzf the nr x nr
// covariance of the processes, Cap = I + XTf*Vzf*t(XTf) (cholCap = chol(Cap), lower) and P = inv(Cap):
//   W = P*X,  Schur(A1) = t(X)*W + inv(Vbeta),  G = Vzf*t(XTf) (nr x n),
//   D1inv = inv(t(XTf)*XTf + inv(Vzf)) = Vzf - G*P*t(G)   (Woodbury, through Cap),
//   D1invB1 = D1inv*t(XTf)*X = G*W,  DinvB_np = W*inv(Schur(A1)),
//   DinvB_nrn = D1inv*(t(XTf) - t(XTf)*X*t(DinvB_np)) = G*P - D1invB1*t(DinvB_np),
//   Schur(A) = I/sigmaSqxi + Q,  Q = P - W*inv(Schur(A1))*t(W).
// Outputs as before: D1inv (nr x nr, full), D1invB1 (nr x p), cholSchurA1_pp = chol(Schur(A1)), DinvB_np (n x p),
// DinvB_nrn (nr x n), cholSchurA_nn = chol(Schur(A)) (lower triangles). tmp_nnr is nr x n workspace.
// XtX and XTildetX are not needed by this formulation.
void primingGLMvc(int n, int p, int r, double *X, double *XTilde, double *XtX, double *XTildetX,
                  double *VBetaInv, double *Vz, std::string &processtype, double *cholCap, double sigmaSqxi,
                  double *tmp_nnr, double *D1inv, double *D1invB1, double *cholSchurA1_pp,
                  double *DinvB_np, double *DinvB_nrn, double *cholSchurA_nn){

  int np = n * p;
  int pp = p * p;
  int nn = n * n;
  int nr = n * r;
  int nnr = n * nr;
  int i = 0, j = 0, k = 0;

  int info = 0;
  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  char const *lside = "L";
  char const *rside = "R";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;
  const double sigmaSqxiInv = 1.0 / sigmaSqxi;

  // W = inv(Cap)*X, held in DinvB_np until it is overwritten below
  F77_NAME(dcopy)(&np, X, &incOne, DinvB_np, &incOne);                                                                  // DinvB_np = X
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholCap, &n, DinvB_np, &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholCap, &n, DinvB_np, &n FCONE FCONE FCONE FCONE);         // DinvB_np = W = inv(Cap)*X

  // chol(Schur(A1)), Schur(A1) = t(X)*W + inv(Vbeta)
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, DinvB_np, &n, &zero, cholSchurA1_pp, &p FCONE FCONE);         // t(X)*W
  F77_NAME(daxpy)(&pp, &one, VBetaInv, &incOne, cholSchurA1_pp, &incOne);                                              // SchurA1 = t(X)*W + VBetaInv
  F77_NAME(dpotrf)(lower, &p, cholSchurA1_pp, &p, &info FCONE); if(info != 0){perror("c++ error: cholSchurA1_pp dpotrf failed\n");}   // chol(Schur(A1))

  // G = Vzf*t(XTf) (nr x n), exploiting the sparsity of XTf; rmul_Vz_XTildeT accumulates, so start from zero
  zeros(tmp_nnr, nnr);
  rmul_Vz_XTildeT(n, r, XTilde, Vz, tmp_nnr, processtype);                                                              // tmp_nnr = G

  // D1invB1 = G*W
  F77_NAME(dgemm)(ntran, ntran, &nr, &p, &n, &one, tmp_nnr, &nr, DinvB_np, &n, &zero, D1invB1, &nr FCONE FCONE);       // D1invB1 = G*W

  // M = G*t(Lcapinv), D1inv = Vzf - M*t(M), and then G*P = M*Lcapinv, in DinvB_nrn
  F77_NAME(dcopy)(&nnr, tmp_nnr, &incOne, DinvB_nrn, &incOne);                                                         // DinvB_nrn = G
  F77_NAME(dtrsm)(rside, lower, ytran, nunit, &nr, &n, &one, cholCap, &n, DinvB_nrn, &nr FCONE FCONE FCONE FCONE);     // DinvB_nrn = M = G*t(Lcapinv)
  F77_NAME(dsyrk)(lower, ntran, &nr, &n, &negone, DinvB_nrn, &nr, &zero, D1inv, &nr FCONE FCONE);                      // D1inv = - M*t(M) (lower)
  for(j = 0; j < nr; j++){
    for(i = j + 1; i < nr; i++){
      D1inv[i*nr + j] = D1inv[j*nr + i];                                                                               // full symmetric - G*P*t(G)
    }
  }
  F77_NAME(dtrsm)(rside, lower, ntran, nunit, &nr, &n, &one, cholCap, &n, DinvB_nrn, &nr FCONE FCONE FCONE FCONE);     // DinvB_nrn = M*Lcapinv = G*P

  // add Vzf to D1inv, to get Vzf - Vzf*t(XTilde)*inv(I + XTilde*Vzf*t(XTilde))*XTilde*Vzf
  // which is equal to inv(t(XTilde)*XTilde + inv(Vzf)) by Sherman-Woodbury-Morrison identity
  if(processtype == "independent.shared" || processtype == "multivariate"){
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
          D1inv[i*n*nr + j*nr + (i*n + k)] += Vz[j*n + k];
        }
      }
    }
  }else if(processtype == "independent"){
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
          D1inv[i*n*nr + j*nr + (i*n + k)] += Vz[i*nn + j*n + k];
        }
      }
    }
  }else if(processtype == "multivariate2"){
    int nrnr = nr * nr;
    F77_NAME(daxpy)(&nrnr, &one, Vz, &incOne, D1inv, &incOne);
  }

  // DinvB_np = W*inv(Schur(A1)) via T = W*t(L1inv); Q = P - T*t(T)
  F77_NAME(dtrsm)(rside, lower, ytran, nunit, &n, &p, &one, cholSchurA1_pp, &p, DinvB_np, &n FCONE FCONE FCONE FCONE); // DinvB_np = T = W*t(L1inv)
  F77_NAME(dcopy)(&nn, cholCap, &incOne, cholSchurA_nn, &incOne);
  F77_NAME(dpotri)(lower, &n, cholSchurA_nn, &n, &info FCONE); if(info != 0){perror("c++ error: Cap dpotri failed\n");}  // cholSchurA_nn = P = inv(Cap)
  F77_NAME(dsyrk)(lower, ntran, &n, &p, &negone, DinvB_np, &n, &one, cholSchurA_nn, &n FCONE FCONE);                    // cholSchurA_nn = Q (lower)
  F77_NAME(dtrsm)(rside, lower, ntran, nunit, &n, &p, &one, cholSchurA1_pp, &p, DinvB_np, &n FCONE FCONE FCONE FCONE); // DinvB_np = W*inv(Schur(A1))

  // DinvB_nrn = G*P - D1invB1*t(DinvB_np)
  F77_NAME(dgemm)(ntran, ytran, &nr, &n, &p, &negone, D1invB1, &nr, DinvB_np, &n, &one, DinvB_nrn, &nr FCONE FCONE);

  // chol(Schur(A)), Schur(A) = I/sigmaSqxi + Q
  for(i = 0; i < n; i++){
    cholSchurA_nn[i*n + i] += sigmaSqxiInv;
  }
  F77_NAME(dpotrf)(lower, &n, cholSchurA_nn, &n, &info FCONE); if(info != 0){perror("c++ error: cholSchurA_nn dpotrf failed\n");}   // chol(Schur(A))

}

// Function for triangular-solve of a vector with only one non-zero entry
void dtrsv_sparse1(double *L, double b, double *x, int n, int k){

  int i = 0, j = 0;
  double sum = 0.0;

  zeros(x, n);

  // Solve for x[k] directly
  x[k] = b / L[k * n + k];

  // Forward solve for x[i] (i > k)
  for(i = k + 1; i < n; i++){

    sum = 0.0;

    // Compute the sum L[i,j] * x[j] for j = k to i - 1
    for(j = k; j < i; j++){
      sum += L[j * n + i] * x[j];
    }

    // Compute x[i]
    x[i] = - sum / L[i * n + i];

  }

}

// Function for the projection for GLM in varying-coefficients model
void projGLMvc(int n, int p, int r, double *X, double *XTilde, double sigmaSqxi, double *Lbeta,
               double *cholVz, std::string &processtype, double *v_eta, double *v_xi, double *v_beta, double *v_z,
               double *D1inv, double *D1invB1, double *cholSchurA1_pp,
               double *DinvB_np, double *DinvB_nrn, double *cholSchurA_nn,
               double *tmp_nr){

  int i = 0;
  int nn = n * n;
  int nr = n * r;

  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;
  const double sigmaxiInv = 1.0 / sqrt(sigmaSqxi);

  // Find components of t(H)*v, where (2n+p+nr)x1 vector v = [v_eta, v_xi, v_beta, v_z]
  F77_NAME(dscal)(&n, &sigmaxiInv, v_xi, &incOne);                                                    // v_xi = v_xi/sigmasqxi
  F77_NAME(daxpy)(&n, &one, v_eta, &incOne, v_xi, &incOne);                                           // v_xi = v_eta + v_xi/sigmasqxi

  F77_NAME(dtrsv)(lower, ytran, nunit, &p, Lbeta, &p, v_beta, &incOne FCONE FCONE FCONE);             // v_beta = LbetatInv*v_beta
  F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, v_eta, &incOne, &one, v_beta, &incOne FCONE);           // v_beta = Xt*v_eta + LbetaInv*v_beta

  if(processtype == "independent.shared" || processtype == "multivariate"){
    for(i = 0; i < r; i++){
      F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVz, &n, &v_z[i*n], &incOne FCONE FCONE FCONE);
    }
  }else if(processtype == "independent"){
    for(i = 0; i < r; i++){
      F77_NAME(dtrsv)(lower, ytran, nunit, &n, &cholVz[i*nn], &n, &v_z[i*n], &incOne FCONE FCONE FCONE);
    }
  }else if(processtype == "multivariate2"){
    // can it be made efficient by passing cholR and chol_iwScale, instead of cholVz? [probably not.]
    F77_NAME(dtrsv)(lower, ytran, nunit, &nr, cholVz, &nr, v_z, &incOne FCONE FCONE FCONE);
  }

  lmulv_XTilde_VC(ytran, n, r, XTilde, v_eta, tmp_nr);
  F77_NAME(daxpy)(&nr, &one, tmp_nr, &incOne, v_z, &incOne);                                          // v_z = XTildet*v_eta + t(LzInv)*v_z

  F77_NAME(dgemv)(ytran, &nr, &n, &one, DinvB_nrn, &nr, v_z, &incOne, &zero, tmp_nr, &incOne FCONE);  // tmp_nr = (BtDInv*v2)_1
  F77_NAME(dgemv)(ntran, &n, &p, &one, DinvB_np, &n, v_beta, &incOne, &one, tmp_nr, &incOne FCONE);   // tmp_nr = BtDInv*v2
  F77_NAME(dscal)(&n, &negone, tmp_nr, &incOne);                                                      // tmp_nr = - BtDInv*v2
  F77_NAME(daxpy)(&n, &one, tmp_nr, &incOne, v_xi, &incOne);                                          // v_xi = v_xi - BtDInv*v2
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholSchurA_nn, &n, v_xi, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholSchurA_nn, &n, v_xi, &incOne FCONE FCONE FCONE);       // v_xi = inv(schurA)*(v_xi - BtDInv*v2)

  F77_NAME(dgemv)(ytran, &nr, &p, &negone, D1invB1, &nr, v_z, &incOne, &zero, tmp_nr, &incOne FCONE); // tmp_nr = - t(B1)*inv(D1)*v_z
  F77_NAME(daxpy)(&p, &one, tmp_nr, &incOne, v_beta, &incOne);                                        // v_beta = v_beta - t(B1)*inv(D1)*v_z
  F77_NAME(dtrsv)(lower, ntran, nunit, &p, cholSchurA1_pp, &p, v_beta, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &p, cholSchurA1_pp, &p, v_beta, &incOne FCONE FCONE FCONE);    // v_beta = inv(schurA1)*(v21 - t(B1)*inv(D1)*v22)

  F77_NAME(dgemv)(ntran, &nr, &p, &negone, D1invB1, &nr, v_beta, &incOne, &zero, tmp_nr, &incOne FCONE);  // tmp_nr = - inv(D1)*B1*v_beta
  F77_NAME(dgemv)(ntran, &nr, &nr, &one, D1inv, &nr, v_z, &incOne, &one, tmp_nr, &incOne FCONE);       // tmp_nr = inv(D1)*v_z - inv(D1)*B1*v_beta
  F77_NAME(dcopy)(&nr, tmp_nr, &incOne, v_z, &incOne);                                                 // v_z = inv(D1)*v_z - inv(D1)*B1*v_beta

  // Find inv(D)*v2 - inv(D)*B*inv(schurA1)*(v1 - t(B)*inv(D)*v2)
  F77_NAME(dgemv)(ytran, &n, &p, &negone, DinvB_np, &n, v_xi, &incOne, &one, v_beta, &incOne FCONE);
  F77_NAME(dgemv)(ntran, &nr, &n, &negone, DinvB_nrn, &nr, v_xi, &incOne, &one, v_z, &incOne FCONE);

}

// Find Kronecker product of two matrices A (rxr) and B (nxn), and store in C (nrxnr)
void kronecker(int r, int n, double *A, double *B, double *C){

  int nr = n * r;
  int i = 0, j = 0, k = 0, l = 0;
  double a_ij = 0.0;
  int C_base = 0;                     // Compute block starting position in C

  for(j = 0; j < r; j++){
    for(i = 0; i < r; i++){

      // A(i, j) in column-major order
      a_ij = A[i + j * r];

      // Compute block starting position in C
      C_base = (j * n) * nr + (i * n);

      // Insert scaled B into the corresponding block of C
      for(l = 0; l < n; l++){
        for(k = 0; k < n; k++){
          C[C_base + l * nr + k] = a_ij * B[k + l * n];
        }
      }
    }
  }
}

// Find Cholesky factor of Kronecker product from the individual Cholesky factors
// Cholesky factors are always considered to be lower-triangular (just kronecker of two triangular matrices)
void chol_kron(int r, int n, double *cholA, double *cholB, double *cholC){

  int nr = n * r;
  int i = 0, j = 0, k = 0, l = 0;
  double a_ij = 0.0;
  int C_base = 0; // Compute block starting position in C

  for(j = 0; j < r; j++){
    for(i = j; i < r; i++){

      // A(i, j) in column-major order
      a_ij = cholA[i + j * r];

      // Compute block starting position in C
      C_base = (j * n) * nr + (i * n);

      // Insert scaled B into the corresponding block of C
      for(l = 0; l < n; l++){
        for(k = l; k < n; k++){
          cholC[C_base + l * nr + k] = a_ij * cholB[k + l * n];
        }
      }
    }
  }

}
// In-place rank-1 downdate of a lower Cholesky factor: L <- chol(L*t(L) - v*t(v)), Krause and Igel (2015)
// with alpha = 1, beta = -1. Column j of the result depends only on column j of the input, so the update
// can overwrite L (the old diagonal entry is kept in ljj). Only the lower triangle is read and written.
// w is n x 1 workspace. Returns 0 on success, or j + 1 if L*t(L) - v*t(v) is not numerically positive
// definite at column j; L is then partly overwritten and must not be used.
int cholRankOneDowndate(int n, double *L, double *v, double *w){

  int j, k;
  const int incOne = 1;
  double b = 1.0, gamma = 0.0, ljj = 0.0, ljjsq = 0.0, wjsq = 0.0, newsq = 0.0, njj = 0.0, tmp = 0.0;

  F77_NAME(dcopy)(&n, v, &incOne, w, &incOne);

  for(j = 0; j < n; j++){

    ljj = L[j*n + j];
    ljjsq = ljj * ljj;
    wjsq = w[j] * w[j];
    gamma = ljjsq * b - wjsq;                // gamma = L[jj]^2*b - w[j]^2
    newsq = ljjsq - wjsq / b;                // new L[jj]^2 = gamma/b
    if(!(gamma > 0.0) || !(newsq > 0.0)){    // also catches NaN
      return j + 1;
    }
    njj = sqrt(newsq);
    L[j*n + j] = njj;

    for(k = j + 1; k < n; k++){
      tmp = L[j*n + k] / ljj;                // old L[kj]/L[jj]
      w[k] -= w[j] * tmp;
      L[j*n + k] = njj * (tmp - w[j] * w[k] / gamma);
    }

    b -= wjsq / ljjsq;

  }

  return 0;

}

// Fast pre-processing for projGLM() on the data with the contiguous block J = [del_start, del_end]
// (0-based, k = del_end - del_start + 1 sites; k = 1 gives leave-one-out) deleted, obtained from the
// full-data outputs of cholSchurGLM() by the partitioned inverse identity
//   inv(M[K,K]) = Sigma[K,K] - Sigma[K,J]*inv(Sigma[J,J])*Sigma[J,K],   Sigma = inv(M), K = complement of J,
// applied with M = Vz + I (Sigma = P) and M = Vz + I + X*Vbeta*t(X) (Sigma = Q). With W = P*X:
//   D1invX_{-J}   = D1invX[K,] + P[K,J]*inv(P[J,J])*W[J,],          W_{-J} = X[K,] - D1invX_{-J},
//   S1_{-J}       = t(X[K,])*W_{-J} + inv(Vbeta),                    DinvB_np_{-J} = W_{-J}*inv(S1_{-J}),
//   I/sigmaSqxi + Q_{-J} = (I/sigmaSqxi + Q)[K,K] - U*t(U),  U = Q[K,J]*t(inv(chol(Q[J,J]))),
// where P[,J] is found by two triangular solves with chol(Vz + I), W[J,] = t(P[,J])*X and
// Q[,J] = P[,J] - DinvB_np*t(W[J,]). S1_{-J} and DinvB_np_{-J} are formed directly from W_{-J}, as in
// cholSchurGLM(), rather than by downdating S1 and DinvB_np: the downdated forms cancel when the
// deleted block is large. The last line is applied as k rank-1 downdates of
// cholSchurDel_n = chol((I/sigmaSqxi + Q)[K,K]) (from cholRowDelUpdate/cholBlockDelUpdate of cholSchur_n),
// in place; every intermediate matrix is >= I/sigmaSqxi, so the downdates are well conditioned.
// Inputs:  X (n x p), cholVzPlusI, D1invX, DinvB_np (n x p) of the full data, VbetaInv (p x p).
// Outputs: D1invX_out, DinvB_np_out ((n-k) x p), cholSchur_p_out (p x p), cholSchurDel_n ((n-k) x (n-k)).
// Workspace: PJ, QJ (n x k), tmp_np (n x p), LP, LQ (k x k), Z (k x p), u, w (n x 1).
// Returns 0 on success; a nonzero value means a factorization or downdate was not numerically positive
// definite, and the caller should recompute the outputs directly (cholSchurGLM on the reduced data).
int cholSchurGLMdel(int n, int p, int del_start, int del_end, double *X, double *cholVzPlusI,
                    double *D1invX, double *DinvB_np, double *VbetaInv,
                    double *D1invX_out, double *DinvB_np_out, double *cholSchur_p_out, double *cholSchurDel_n,
                    double *PJ, double *QJ, double *tmp_np, double *LP, double *LQ, double *Z,
                    double *u, double *w){

  int l;
  int info = 0;
  int k = del_end - del_start + 1;
  int nk = n * k;
  int np = n * p;
  int pp = p * p;
  int nMs = n - del_start;
  int nnk = n - k;
  int nnkp = nnk * p;
  char const *lower = "L";
  char const *ytran = "T";
  char const *ntran = "N";
  char const *nunit = "N";
  char const *lside = "L";
  char const *rside = "R";
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;

  // PJ = P[,J] = inv(Vz + I)*E_J; rows above del_start of inv(L)*E_J are zero, so the forward solve
  // starts at row del_start
  zeros(PJ, nk);
  for(l = 0; l < k; l++){
    PJ[l*n + del_start + l] = 1.0;
  }
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &nMs, &k, &one, &cholVzPlusI[del_start*n + del_start], &n, &PJ[del_start], &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &k, &one, cholVzPlusI, &n, PJ, &n FCONE FCONE FCONE FCONE);   // PJ = P[,J]

  // LP = chol(P[J,J])
  copyMatrixRowBlock(PJ, n, k, LP, del_start, del_end);
  F77_NAME(dpotrf)(lower, &k, LP, &k, &info FCONE); if(info != 0){return 1;}

  // Z = W[J,] = t(P[,J])*X (k x p)
  F77_NAME(dgemm)(ytran, ntran, &k, &p, &n, &one, PJ, &n, X, &n, &zero, Z, &k FCONE FCONE);

  // QJ = Q[,J] = P[,J] - DinvB_np*t(W[J,])
  F77_NAME(dcopy)(&nk, PJ, &incOne, QJ, &incOne);
  F77_NAME(dgemm)(ntran, ytran, &n, &k, &p, &negone, DinvB_np, &n, Z, &k, &one, QJ, &n FCONE FCONE);

  // D1invX_out = (D1invX + P[,J]*inv(P[J,J])*W[J,])[K,]
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &k, &p, &one, LP, &k, Z, &k FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &k, &p, &one, LP, &k, Z, &k FCONE FCONE FCONE FCONE);             // Z = inv(P[J,J])*W[J,]
  F77_NAME(dcopy)(&np, D1invX, &incOne, tmp_np, &incOne);
  F77_NAME(dgemm)(ntran, ntran, &n, &p, &k, &one, PJ, &n, Z, &k, &one, tmp_np, &n FCONE FCONE);
  copyMatrixDelRowBlock(tmp_np, n, p, D1invX_out, del_start, del_end);

  // W_{-J} = X[K,] - D1invX_{-J} (in tmp_np, leading dimension n - k); X[K,] held in DinvB_np_out for now
  copyMatrixDelRowBlock(X, n, p, DinvB_np_out, del_start, del_end);                                              // DinvB_np_out = X[K,]
  F77_NAME(dcopy)(&nnkp, DinvB_np_out, &incOne, tmp_np, &incOne);
  F77_NAME(daxpy)(&nnkp, &negone, D1invX_out, &incOne, tmp_np, &incOne);                                         // tmp_np = W_{-J}

  // cholSchur_p_out = chol(t(X[K,])*W_{-J} + inv(Vbeta))
  F77_NAME(dcopy)(&pp, VbetaInv, &incOne, cholSchur_p_out, &incOne);
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &nnk, &one, DinvB_np_out, &nnk, tmp_np, &nnk, &one, cholSchur_p_out, &p FCONE FCONE);
  F77_NAME(dpotrf)(lower, &p, cholSchur_p_out, &p, &info FCONE); if(info != 0){return 2;}

  // DinvB_np_out = W_{-J}*inv(S1_{-J})
  F77_NAME(dcopy)(&nnkp, tmp_np, &incOne, DinvB_np_out, &incOne);
  F77_NAME(dtrsm)(rside, lower, ytran, nunit, &nnk, &p, &one, cholSchur_p_out, &p, DinvB_np_out, &nnk FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(rside, lower, ntran, nunit, &nnk, &p, &one, cholSchur_p_out, &p, DinvB_np_out, &nnk FCONE FCONE FCONE FCONE);

  // chol(I/sigmaSqxi + Q_{-J}): k rank-1 downdates by the columns of U = Q[K,J]*t(inv(chol(Q[J,J])))
  copyMatrixRowBlock(QJ, n, k, LQ, del_start, del_end);
  F77_NAME(dpotrf)(lower, &k, LQ, &k, &info FCONE); if(info != 0){return 3;}                                     // LQ = chol(Q[J,J])
  F77_NAME(dtrsm)(rside, lower, ytran, nunit, &n, &k, &one, LQ, &k, QJ, &n FCONE FCONE FCONE FCONE);            // QJ = Q[,J]*t(inv(LQ))
  for(l = 0; l < k; l++){
    copyVecExcludingBlock(&QJ[l*n], u, n, del_start, del_end);                                                   // u = U[,l]
    info = cholRankOneDowndate(nnk, cholSchurDel_n, u, w);
    if(info != 0){return 4;}
  }

  return 0;

}
