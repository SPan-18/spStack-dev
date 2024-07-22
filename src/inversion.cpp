#define USE_FC_LEN_T
#include <string>
#include "inversion.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

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

  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, v2, &incOne, &zero, tmp_n1, &incOne FCONE);  // tmp_n1 = Vz*v2
  F77_NAME(dcopy)(&n, tmp_n1, &incOne, tmp_n2, &incOne);  // tmp_n1 = tmp_n2
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n2, &incOne FCONE FCONE FCONE);  // tmp_n2 = VyInv*Vz*v2
  F77_NAME(dgemv)(ntran, &n, &n, &negone, Vz, &n, tmp_n2, &incOne, &one, tmp_n1, &incOne FCONE);  // tmp_n1 = inv(VzInv+deltasq*I)*v2
  F77_NAME(dcopy)(&n, tmp_n1, &incOne, out_n, &incOne);  // out_n = tmp_n1 = inv(D)*v2
  F77_NAME(dcopy)(&p, v1, &incOne, tmp_p1, &incOne);  // tmp_p1 = v1
  F77_NAME(dgemv)(ytran, &n, &p, &negdeltasqInv, X, &n, tmp_n1, &incOne, &one, tmp_p1, &incOne FCONE);  // tmp_p1 = v1 - t(B)*inv(D)*v2

  F77_NAME(dcopy)(&pp, VbetaInv, &incOne, tmp_pp, &incOne);  // tmp_pp = VbetaInv
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &deltasqInv, X, &n, X, &n, &one, tmp_pp, &p FCONE FCONE);  // tmp_pp = A = (1/deltasq)*XtX+VbetaInv

  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &deltasqInv, Vz, &n, X, &n, &zero, tmp_np1, &n FCONE FCONE);  // tmp_np1 = Vz*B
  F77_NAME(dcopy)(&np, tmp_np1, &incOne, tmp_np2, &incOne);  // tmp_np2 = tmp_np1 = Vz*B
  F77_NAME(dtrsm)(lside, lower, ntran, nunit, &n, &p, &one, cholVy, &n, tmp_np2, &n FCONE FCONE FCONE FCONE); // tmp_np2 = LyInv*Vz*B
  F77_NAME(dtrsm)(lside, lower, ytran, nunit, &n, &p, &one, cholVy, &n, tmp_np2, &n FCONE FCONE FCONE FCONE); // tmp_np2 = VyInv*Vz*B
  F77_NAME(dgemm)(ntran, ntran, &n, &p, &n, &negone, Vz, &n, tmp_np2, &n, &one, tmp_np1, &n FCONE FCONE);  // tmp_np1 = (Vz - VzVyinv*Vz)*B = inv(D)*B
  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &negdeltasqInv, X, &n, tmp_np1, &n, &one, tmp_pp, &p FCONE FCONE);  // tmp_pp = Schur(A) = A - t(B)*inv(D)*B
  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}  // chol(Schur(A))
  F77_NAME(dtrsv)(lower, ntran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);  // tmp_p1 = inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dcopy)(&p, tmp_p1, &incOne, out_p, &incOne);  // out_p = first p elements of Mv

  F77_NAME(dgemv)(ntran, &n, &p, &deltasqInv, X, &n, tmp_p1, &incOne, &zero, tmp_n1, &incOne FCONE);  // tmp_n1 = B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dgemv)(ntran, &n, &n, &one, Vz, &n, tmp_n1, &incOne, &zero, tmp_n2, &incOne FCONE);  // tmp_n2 = Vz * tmp_n1
  F77_NAME(dcopy)(&n, tmp_n2, &incOne, tmp_n1, &incOne);  // tmp_n1 = tmp_n2
  F77_NAME(dtrsv)(lower, ntran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);
  F77_NAME(dtrsv)(lower, ytran, nunit, &n, cholVy, &n, tmp_n1, &incOne FCONE FCONE FCONE);  // tmp_n1 = VyInv*Vz*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(dgemv)(ntran, &n, &n, &negone, Vz, &n, tmp_n1, &incOne, &one, tmp_n2, &incOne FCONE);  // tmp_n2 = inv(D)*B*inv(Schur(A))*(v1-BtDinvB)
  F77_NAME(daxpy)(&n, &negone, tmp_n2, &incOne, out_n, &incOne); // out_n = inv(D)*v2 - inv(D)*B*inv(Schur(A))*(v1-BtDinvB)

}
