#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "MatrixAlgos.h"

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

// Cholesky factor update after deletion of a row/column
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
    w = (double *)R_chk_realloc(w, nk * sizeof(double));
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

    }  // End rank-one update for first row/column deletion

  }else if(0 < del < n - 1){

    int delPlusOne = del + 1;
    int indexL1 = 0;

    nk = n - delPlusOne;

    copySubmat(L, n, n, L1, n1, n1, 0, 0, 0, 0, del, del);
    copySubmat(L, n, n, L1, n1, n1, delPlusOne, 0, del, 0, nk, del);

    w = (double *)R_chk_realloc(w, nk * sizeof(double));
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


  }else{
    perror("Row/column deletion index out of bounds.");
  }

}

// No memory allocation inside 'hot' function
void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1,
                 double *tmp_pp, double *tmp_np1, double*tmp_np2,
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
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}    // chol(Schur(A))
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

  double *tmp_n1 = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n1, n);
  double *tmp_n2 = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n2, n);

  double *tmp_np1 = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np1, np);
  double *tmp_np2 = (double *) R_alloc(np, sizeof(double)); zeros(tmp_np2, np);

  double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);

  double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);

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
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: dpotrf failed\n");}    // chol(Schur(A))
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
