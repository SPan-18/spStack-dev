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
  if(del_start > del_end || del_start == del_end){
    perror("Block Start index must be at least 1 less than End index.");
  }
  if(del_start < 0 || del_end > n){
    perror("Block index to delete is out of bounds.");
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
void cholSchurGLM(double *X, int n, int p, double sigmaSqxi, double *XtX, double *VbetaInv,
                  double *Vz, double *cholVzPlusI, double *tmp_nn, double *tmp_np,
                  double *tmp_pn, double *tmp_nn2, double *out_pp, double *out_nn, double *D1invB1){

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
  const double one = 1.0;
  const double negone = -1.0;
  const double zero = 0.0;
  const int incOne = 1;
  const double sigmaSqxi2 = (sigmaSqxi + 1.0) / sigmaSqxi;

  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &one, Vz, &n, X, &n, &zero, tmp_np, &n FCONE FCONE);                  // tmp_np = Vz*X
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholVzPlusI, &n, tmp_np, &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholVzPlusI, &n, tmp_np, &n FCONE FCONE FCONE FCONE); // tmp_np = inv(VzInv+I)*X
  F77_NAME(dcopy)(&np, tmp_np, &incOne, D1invB1, &incOne);
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np, &n, &zero, out_pp, &p FCONE FCONE);              // out_pp = t(X)*inv(VzInv+I)*X
  F77_NAME(dscal)(&pp, &negone, out_pp, &incOne);                                                                 // out_pp = - XtD1invX
  F77_NAME(daxpy)(&pp, &one, XtX, &incOne, out_pp, &incOne);                                                      // out_pp = XtX - XtD1invX
  F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, out_pp, &incOne);                                                 // out_pp = XtX + VbetaInv - XtD1invX
  F77_NAME(dpotrf)(lower, &p, out_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}      // chol(Schur(A1))
  F77_NAME(daxpy)(&np, &negone, X, &incOne, tmp_np, &incOne);                                                     // tmp_np = inv(VzInv+I)*X-X
  F77_NAME(dscal)(&np, &negone, tmp_np, &incOne);                                                                 // tmp_np = X-inv(VzInv+I)*X
  transpose_matrix(tmp_np, n, p, tmp_pn);                                                                         // tmp_pn = t(tmp_np)
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &p, &n, &one, out_pp, &p, tmp_pn, &p FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &p, &n, &one, out_pp, &p, tmp_pn, &p FCONE FCONE FCONE FCONE);      // tmp_pn = inv(schur(A1))*(X-inv(VzInv+I)*X) and RETURN
  F77_NAME(dgemm)(ntran, ntran, &n, &n, &p, &one, X, &n, tmp_pn, &p, &zero, tmp_nn, &n FCONE FCONE);              // tmp_nn = X*inv(schur(A1))*(X-inv(VzInv+I)*X)
  F77_NAME(dscal)(&nn, &negone, tmp_nn, &incOne);                                                                 // tmp_nn = -X*inv(schur(A1))*(X-inv(VzInv+I)*X)

  for(i = 0; i < n; i++){
    tmp_nn[i * n + i] += 1.0;
  }

  F77_NAME(dgemm)(ntran, ntran, &n, &n, &n, &one, Vz, &n, tmp_nn, &n, &zero, out_nn, &n FCONE FCONE);             // out_nn = Vz*(I-X*inv(schur(A1))*(X-inv(VzInv+I)*X))
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &n, &one, cholVzPlusI, &n, out_nn, &n FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &n, &one, cholVzPlusI, &n, out_nn, &n FCONE FCONE FCONE FCONE); // out_nn = inv(D1)*out_nn
  F77_NAME(dcopy)(&nn, out_nn, &incOne, tmp_nn2, &incOne);                                                        // return DinvB_nn

  for(i = 0; i < n; i++){
    tmp_nn[i * n + i] -= 1.0;                                                                                     // tmp_nn = -X*inv(schur(A1))*(X-inv(VzInv+I)*X)
  }

  F77_NAME(dscal)(&nn, &negone, tmp_nn, &incOne);                                                                 // tmp_nn = X*inv(schur(A1))*(X-inv(VzInv+I)*X)
  F77_NAME(daxpy)(&nn, &one, tmp_nn, &incOne, out_nn, &incOne);                                                   // out_nn = X*tmp_pn + out_nn
  F77_NAME(dscal)(&nn, &negone, out_nn, &incOne);                                                                 // out_nn = - BtDinvB

  for(i = 0; i < n; i++){
    out_nn[i * n + i] += sigmaSqxi2;                                                                              // out_nn = A - BtDinvB
  }

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

// memory allocation inside hot function
void inversionLM2(double *X, int n, int p, double deltasq, double *VbetaInv,
                  double *Vz, double *cholVy, double *v1, double *v2,
                  double *out_p, double *out_n){

  int pp = p * p;
  int np = n * p;

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

  double *tmp_n1 = (double *) R_chk_calloc(n, sizeof(double)); zeros(tmp_n1, n);
  double *tmp_n2 = (double *) R_chk_calloc(n, sizeof(double)); zeros(tmp_n2, n);

  double *tmp_np1 = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np1, np);
  double *tmp_np2 = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np2, np);

  double *tmp_p1 = (double *) R_chk_calloc(p, sizeof(double)); zeros(tmp_p1, p);

  double *tmp_pp = (double *) R_chk_calloc(pp, sizeof(double)); zeros(tmp_pp, pp);

  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, v2, &incOne, &zero, tmp_n1, &incOne FCONE);                     // tmp_n1 = Vz*v2
  F77_NAME(dcopy)(&n, tmp_n1, &incOne, tmp_n2, &incOne);                                                       // tmp_n1 = tmp_n2
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);                     // tmp_n2 = VyInv*Vz*v2
  F77_NAME(dgemv)(ntran, &n, &n, &negone, Vz, &n, tmp_n2, &incOne, &one, tmp_n1, &incOne FCONE);               // tmp_n1 = inv(VzInv+deltasq*I)*v2
  F77_NAME(dcopy)(&n, tmp_n1, &incOne, out_n, &incOne);                                                        // out_n = tmp_n1 = inv(D)*v2
  F77_NAME(dcopy)(&p, v1, &incOne, tmp_p1, &incOne);                                                           // tmp_p1 = v1
  F77_NAME(dgemv)(ytran, &n, &p, &negdeltasqInv, X, &n, tmp_n1, &incOne, &one, tmp_p1, &incOne FCONE);         // tmp_p1 = v1 - t(B)*inv(D)*v2

  F77_NAME(dcopy)(&pp, VbetaInv, &incOne, tmp_pp, &incOne);                                                    // tmp_pp = VbetaInv
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &deltasqInv, X, &n, X, &n, &one, tmp_pp, &p FCONE FCONE);          // tmp_pp = A = (1/deltasq)*XtX+VbetaInv

  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &deltasqInv, Vz, &n, X, &n, &zero, tmp_np1, &n FCONE FCONE);       // tmp_np1 = Vz*B
  F77_NAME(dcopy)(&np, tmp_np1, &incOne, tmp_np2, &incOne);                                                    // tmp_np2 = tmp_np1 = Vz*B
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholVy, &n, tmp_np2, &n FCONE FCONE FCONE FCONE);  // tmp_np2 = LyInv*Vz*B
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholVy, &n, tmp_np2, &n FCONE FCONE FCONE FCONE);  // tmp_np2 = VyInv*Vz*B
  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &negone, Vz, &n, tmp_np2, &n, &one, tmp_np1, &n FCONE FCONE);      // tmp_np1 = (Vz - VzVyinv*Vz)*B = inv(D)*B
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &negdeltasqInv, X, &n, tmp_np1, &n, &one, tmp_pp, &p FCONE FCONE); // tmp_pp = Schur(A) = A - t(B)*inv(D)*B
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}   // chol(Schur(A))
  F77_NAME(dtrsv)(lower, ntran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                     // tmp_p1 = inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dcopy)(&p, tmp_p1, &incOne, out_p, &incOne);                                                        // out_p = first p elements of Mv

  F77_NAME(dgemv)(ntran, &n, &p, &deltasqInv, X, &n, tmp_p1, &incOne, &zero, tmp_n1, &incOne FCONE);           // tmp_n1 = B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, tmp_n1, &incOne, &zero, tmp_n2, &incOne FCONE);                 // tmp_n2 = Vz * tmp_n1
  F77_NAME(dcopy)(&n, tmp_n2, &incOne, tmp_n1, &incOne);                                                       // tmp_n1 = tmp_n2
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);                     // tmp_n1 = VyInv*Vz*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dgemv)(ntran, &n, &n, &negone, Vz, &n, tmp_n1, &incOne, &one, tmp_n2, &incOne FCONE);               // tmp_n2 = inv(D)*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(daxpy)(&n, &negone, tmp_n2, &incOne, out_n, &incOne);                                               // out_n = inv(D)*v2 - inv(D)*B*inv(Schur(A))*(v1-BtDinvB)

  R_chk_free(tmp_n1);
  R_chk_free(tmp_n2);
  R_chk_free(tmp_np1);
  R_chk_free(tmp_np2);
  R_chk_free(tmp_p1);
  R_chk_free(tmp_pp);

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
void projGLM(double *X, int n, int p, double *v_eta, double *v_xi, double *v_beta, double *v_z,
             double *cholpSchur, double *cholnSchur, double sigmaSqxi, double *Lbeta, double *Lz,
             double *Vz, double *cholVzPlusI, double *D1invB1, double *DinvBpn, double *DinvBnn,
             double *tmp_n, double *tmp_p){

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

  // Find inv(VzInv+I)*v22
  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, v_z, &incOne, &zero, tmp_n, &incOne FCONE);         // tmp_n = Vz*v22
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVzPlusI, &n, tmp_n, &incOne FCONE FCONE FCONE);     // tmp_n = D1inv*v22

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

  // Find DInvB*inv(schurA1)*(v1 - BtDInvv2)
  F77_NAME(dgemv)(ntran, &p, &n, &one, DinvBpn, &p, v_xi, &incOne, &zero, tmp_p, &incOne FCONE);
  F77_NAME(dgemv)(ntran, &n, &n, &one, DinvBnn, &n, v_xi, &incOne, &zero, tmp_n, &incOne FCONE);   // (tmp_p,n) = (DinvB) * inv(sA)*(v1-BtDInvv2)

  // Find v_beta, v_z
  F77_NAME(daxpy)(&p, &negone, tmp_p, &incOne, v_beta, &incOne);
  F77_NAME(daxpy)(&n, &negone, tmp_n, &incOne, v_z, &incOne);

}

// Function to transpose a matrix in column-major form
void transpose_matrix(double *M, int nrow, int ncol, double *Mt){

  int i, j;

  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      Mt[j + i * ncol] = M[i + j * nrow];
    }
  }

}

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
void primingGLMvc(int n, int p, int r, double *X, double *XTilde, double *XtX, double *XTildetX,
                  double *VBetaInv, double *Vz, int sharedP, double *cholCap, double sigmaSqxi,
                  double *tmp_nnr, double *D1inv, double *D1invB1, double *cholSchurA1_pp,
                  double *DinvB_pn, double *DinvB_nrn, double *cholSchurA_nn){

  int np = n * p;
  int pp = p * p;
  int nn = n * n;
  int nr = n * r;
  int nrp = nr * p;
  int nnr = n * nr;
  int nrnr = nr * nr;
  int i = 0, j = 0, k = 0;

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
  const double sigmaSqxi2 = (sigmaSqxi + 1.0) / sigmaSqxi;

  // F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &one, XtVzXtt, &n, X, &n, &zero, tmp_np1, &n FCONE FCONE);                 // tmp_np1 = X_tilde*Vz*t(X_tilde)*X
  // F77_NAME(dcopy)(&np, tmp_np1, &incOne, tmp_np2, &incOne);                                                            // tmp_n2 = X_tilde*Vz*t(X_tilde)*X
  // F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholCap, &n, tmp_np1, &n FCONE FCONE FCONE FCONE);
  // F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholCap, &n, tmp_np1, &n FCONE FCONE FCONE FCONE);         // tmp_np1 = inv(I+XTildeVzXTildet)*X_tilde*Vz*t(X_tilde)*X
  // F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &negone, XtVzXtt, &n, tmp_np1, &n, &one, tmp_np2, &n FCONE FCONE);         // tmp_np2 = X_tilde*inv(D1)*t(X_tilde)*X
  // F77_NAME(dscal)(&np, &negone, tmp_np2, &incOne);                                                                     // tmp_np2 = -X_tilde*inv(D1)*t(X_tilde)*X
  // F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np2, &n, &zero, tmp_pp, &p FCONE FCONE);                  // tmp_pp = -t(X)*X_tilde*inv(D1)*t(X_tilde)*X
  // F77_NAME(daxpy)(&pp, &one, XtX, &incOne, tmp_pp, &incOne);                                                           // tmp_pp = XtX - t(X)*X_tilde*inv(D1)*t(X_tilde)*X
  // F77_NAME(daxpy)(&pp, &one, VBetaInv, &incOne, tmp_pp, &incOne);                                                      // tmp_pp = XtX + VbetaInv - t(X)*X_tilde*inv(D1)*t(X_tilde)*X = Schur(A1)
  // F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}           // tmp_pp = chol(Schur(A1))
  // F77_NAME(daxpy)(&np, &one, X, &incOne, tmp_np2, &incOne);                                                            // tmp_np2 = X - X_tilde*inv(D1)*t(X_tilde)*X
  // transpose_matrix(tmp_np2, n, p, tmp_pn);                                                                             // tmp_pn = t(tmp_np2) = t(X - X_tilde*inv(D1)*t(X_tilde)*X)
  // F77_NAME(dtrsm)(lside, lower, ntran, nunit, &p, &n, &one, tmp_pp, &p, tmp_pn, &p FCONE FCONE FCONE FCONE);
  // F77_NAME(dtrsm)(lside, lower, ytran, nunit, &p, &n, &one, tmp_pp, &p, tmp_pn, &p FCONE FCONE FCONE FCONE);           // tmp_np1 = inv(schur(A1))*t(X - X_tilde*inv(D1)*t(X_tilde)*X) and RETURN

  // printMtrx(tmp_pn, p, n);

  // Find inv(D1)
  zeros(D1inv, nrnr);
  // Step 1: Find invchol(I+XTVzXTt)*XT
  for(i = 0; i < r; i++){
    for(j = 0; j < n; j++){
      dtrsv_sparse1(cholCap, XTilde[i*n + j], &D1inv[(i*n + j)*n], n, j);
    }
  }
  if(sharedP){
    for(i = 0; i < r; i++){
      F77_NAME(dgemm)(ntran, ntran, &n, &n, &n, &one, &D1inv[i*nn], &n, Vz, &n, &zero, &tmp_nnr[i*nn], &n FCONE FCONE);
    }
  }else{
    for(i = 0; i < r; i++){
      F77_NAME(dgemm)(ntran, ntran, &n, &n, &n, &one, &D1inv[i*nn], &n, &Vz[i*nn], &n, &zero, &tmp_nnr[i*nn], &n FCONE FCONE);
    }
  }
  F77_NAME(dgemm)(ytran, ntran, &nr, &nr, &n, &one, tmp_nnr, &n, tmp_nnr, &n, &zero, D1inv, &nr FCONE FCONE);
  F77_NAME(dscal)(&nrnr, &negone, D1inv, &incOne);
  if(sharedP){
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
          D1inv[i*n*nr + j*nr + (i*n + k)] += Vz[j*n + k];
        }
      }
    }
  }else{
    for(i = 0; i < r; i++){
      for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
          D1inv[i*n*nr + j*nr + (i*n + k)] += Vz[i*nn + j*n + k];
        }
      }
    }
  }

  F77_NAME(dgemm)(ntran, ntran, &nr, &p, &nr, &one, D1inv, &nr, XTildetX, &nr, &zero, tmp_nnr, &nr FCONE FCONE);       // inv(D1)*XTildetX
  F77_NAME(dcopy)(&nrp, tmp_nnr, &incOne, D1invB1, &incOne);                                                           // D1invB1 = inv(D1)*XTildetX
  lmulm_XTilde_VC(ntran, n, r, p, XTilde, tmp_nnr, DinvB_pn);                                                          // XTilde*inv(D1)*XTildetX
  F77_NAME(dscal)(&np, &negone, DinvB_pn, &incOne);                                                                    // - XTilde*inv(D1)*XTildetX

  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, DinvB_pn, &n, &zero, cholSchurA1_pp, &p FCONE FCONE);         // - t(X)*XTilde*inv(D1)*XTildetX
  F77_NAME(daxpy)(&pp, &one, XtX, &incOne, cholSchurA1_pp, &incOne);                                                   // XtX - t(X)*XTilde*inv(D1)*XTildetX
  F77_NAME(daxpy)(&pp, &one, VBetaInv, &incOne, cholSchurA1_pp, &incOne);                                              // SchurA1 = XtX + VBetaInv - t(X)*XTilde*inv(D1)*XTildetX
  F77_NAME(dpotrf)(lower, &p, cholSchurA1_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}   // chol(Schur(A1))

  F77_NAME(daxpy)(&np, &one, X, &incOne, DinvB_pn, &incOne);                                                           // X - XTilde*inv(D1)*XTildetX
  F77_NAME(dcopy)(&np, DinvB_pn, &incOne, tmp_nnr, &incOne);
  transpose_matrix(tmp_nnr, n, p, DinvB_pn);
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &p, &n, &one, cholSchurA1_pp, &p, DinvB_pn, &p FCONE FCONE FCONE FCONE);
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &p, &n, &one, cholSchurA1_pp, &p, DinvB_pn, &p FCONE FCONE FCONE FCONE); // inv(schurA1)*t(X - XTilde*inv(D1)*XTildetX)

  F77_NAME(dgemm)(ntran, ntran, &nr, &n, &p, &one, XTildetX, &nr, DinvB_pn, &p, &zero, tmp_nnr, &nr FCONE FCONE);      // t(XTilde)*X*DinvB_pn
  F77_NAME(dscal)(&nnr, &negone, tmp_nnr, &incOne);                                                                    // - t(XTilde)*X*DinvB_pn
  addXTildeTransposeToMatrixByRow(XTilde, tmp_nnr, n, r);                                                              // t(XTilde) - t(XTilde)*X*DinvB_pn
  F77_NAME(dgemm)(ntran, ntran, &nr, &n, &nr, &one, D1inv, &nr, tmp_nnr, &nr, &zero, DinvB_nrn, &nr FCONE FCONE);      // inv(D1)*(t(XTilde) - t(XTilde)*X*DinvB_pn)

  lmulm_XTilde_VC(ntran, n, r, n, XTilde, DinvB_nrn, cholSchurA_nn);                                                   // XTilde*DinvB_nrn
  F77_NAME(dgemm)(ntran, ntran, &n, &n, &p, &one, X, &n, DinvB_pn, &p, &one, cholSchurA_nn, &n FCONE FCONE);           // XTilde*DinvB_nrn + X*DinvB_pn
  F77_NAME(dscal)(&nn, &negone, cholSchurA_nn, &incOne);
  for(i = 0; i < n; i++){
    cholSchurA_nn[i*n + i] += sigmaSqxi2;
  }
  F77_NAME(dpotrf)(lower, &n, cholSchurA_nn, &n, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}   // chol(Schur(A))

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
               double *cholVz, int sharedP, double *v_eta, double *v_xi, double *v_beta, double *v_z,
               double *D1inv, double *D1invB1, double *cholSchurA1_pp,
               double *DinvB_pn, double *DinvB_nrn, double *cholSchurA_nn,
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

  if(sharedP){
    for(i = 0; i < r; i++){
      F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVz, &n, &v_z[i*n], &incOne FCONE FCONE FCONE);
    }
  }else{
    for(i = 0; i < r; i++){
      F77_NAME(dtrsv)(lower, ytran, nunit, &n, &cholVz[i*nn], &n, &v_z[i*n], &incOne FCONE FCONE FCONE);
    }
  }
  lmulv_XTilde_VC(ytran, n, r, XTilde, v_eta, tmp_nr);
  F77_NAME(daxpy)(&nr, &one, tmp_nr, &incOne, v_z, &incOne);                                          // v_z = XTildet*v_eta + t(LzInv)*v_z

  F77_NAME(dgemv)(ytran, &nr, &n, &one, DinvB_nrn, &nr, v_z, &incOne, &zero, tmp_nr, &incOne FCONE);  // tmp_nr = (BtDInv*v2)_1
  F77_NAME(dgemv)(ytran, &p, &n, &one, DinvB_pn, &p, v_beta, &incOne, &one, tmp_nr, &incOne FCONE);   // tmp_nr = BtDInv*v2
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
  F77_NAME(dgemv)(ntran, &p, &n, &negone, DinvB_pn, &p, v_xi, &incOne, &one, v_beta, &incOne FCONE);
  F77_NAME(dgemv)(ntran, &nr, &n, &negone, DinvB_nrn, &nr, v_xi, &incOne, &one, v_z, &incOne FCONE);

}