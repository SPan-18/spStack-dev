#include <string>
#include <Rinternals.h>

void cholRankOneUpdate(int n, double *L1, double alpha, double beta,
                       double *v, double *L2, double *w);

void cholRowDelUpdate(int n, double *L, int del, double *L1, double *w);

void cholBlockDelUpdate(int n, double *L, int del_start, int del_end, double *L1, double *tmpL1, double *w);

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

void upperTri_lowerTri(double *M, int n);

void primingGLMvc(int n, int p, int r, double *X, double *XTilde, double *XtX, double *XTildetX,
                  double *VBetaInv, double *Vz, std::string &processtype, double *cholCap, double sigmaSqxi,
                  double *tmp_nnr, double *D1inv, double *D1invB1, double *cholSchurA1_pp,
                  double *DinvB_pn, double *DinvB_nrn, double *cholSchurA_nn);

void dtrsv_sparse1(double *L, double b, double *x, int n, int k);

void projGLMvc(int n, int p, int r, double *X, double *XTilde, double sigmaSqxi, double *Lbeta,
               double *cholVz, std::string &processtype, double *v_eta, double *v_xi, double *v_beta, double *v_z,
               double *D1inv, double *D1invB1, double *cholSchurA1_pp,
               double *DinvB_pn, double *DinvB_nrn, double *cholSchurA_nn,
               double *tmp_nr);

void kronecker(int r, int n, double *A, double *B, double *C);

void chol_kron(int r, int n, double *cholA, double *cholB, double *cholC);