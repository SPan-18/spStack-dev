#include <string>
#include <Rinternals.h>

void cholRankOneUpdate(int n, double *L1, double alpha, double beta,
                       double *v, double *L2, double *w);

void cholRowDelUpdate(int n, double *L, int del, double *L1, double *w);

void cholSchurGLM(double *X, int n, int p, double sigmaSqxi, double *XtX, double *VbetaInv,
                  double *Vz, double *cholVzPlusI, double *tmp_nn, double *tmp_np,
                  double *tmp_pn, double *tmp_nn2, double *out_pp, double *out_nn, double *D1invB1);

void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1,
                 double *tmp_pp, double *tmp_np1,
                 double *outp, double *outn, int LOO);

void inversionLM2(double *X, int n, int p, double deltasq, double *VbetaInv,
                  double *Vz, double *cholVy, double *v1, double *v2,
                  double *out_p, double *out_n);

int mapIndex(int i, int j, int nRowB, int nColB, int startRowB, int startColB, int nRowA);

void projGLM(double *X, int n, int p, double *v_eta, double *v_xi, double *v_beta, double *v_z,
             double *cholpSchur, double *cholnSchur, double sigmaSqxi, double *Lbeta, double *Lz,
             double *Vz, double *cholVzPlusI, double *D1invB1, double *DinvBpn, double *DinvBnn,
             double *tmp_n, double *tmp_p);

void transpose_matrix(double *M, int nrow, int ncol, double *Mt);