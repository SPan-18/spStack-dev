#include <string>
#include <Rinternals.h>

void cholRankOneUpdate(int n, double *L1, double alpha, double beta,
                       double *v, double *L2, double *w);

void cholRowDelUpdate(int n, double *L, int del, double *L1, double *w);

void cholBlockDelUpdate(int n, double *L, int del_start, int del_end, double *L1, double *tmpL1, double *w);

void cholSchurGLM(double *X, int n, int p, double sigmaSqxi, double *VbetaInv, double *cholVzPlusI,
                  double *tmp_np, double *DinvB_np, double *out_pp, double *out_nn, double *D1invB1);

void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1,
                 double *tmp_pp, double *tmp_np1,
                 double *outp, double *outn, int LOO);

int mapIndex(int i, int j, int nRowB, int nColB, int startRowB, int startColB, int nRowA);

void projGLM(double *X, int n, int p, double *v_eta, double *v_xi, double *v_beta, double *v_z,
             double *cholpSchur, double *cholnSchur, double sigmaSqxi, double *Lbeta, double *Lz,
             double *cholVzPlusI, double *D1invB1, double *DinvBnp, double *tmp_n, double *tmp_p);


void upperTri_lowerTri(double *M, int n);

void primingGLMvc(int n, int p, int r, double *X, double *XTilde, double *XtX, double *XTildetX,
                  double *VBetaInv, double *Vz, std::string &processtype, double *cholCap, double sigmaSqxi,
                  double *tmp_nnr, double *D1inv, double *D1invB1, double *cholSchurA1_pp,
                  double *DinvB_np, double *DinvB_nrn, double *cholSchurA_nn);

void dtrsv_sparse1(double *L, double b, double *x, int n, int k);

void projGLMvc(int n, int p, int r, double *X, double *XTilde, double sigmaSqxi, double *Lbeta,
               double *cholVz, std::string &processtype, double *v_eta, double *v_xi, double *v_beta, double *v_z,
               double *D1inv, double *D1invB1, double *cholSchurA1_pp,
               double *DinvB_np, double *DinvB_nrn, double *cholSchurA_nn,
               double *tmp_nr);

int cholRankOneDowndate(int n, double *L, double *v, double *w);

int cholSchurGLMdel(int n, int p, int del_start, int del_end, double *X, double *cholVzPlusI,
                    double *D1invX, double *DinvB_np, double *VbetaInv,
                    double *D1invX_out, double *DinvB_np_out, double *cholSchur_p_out, double *cholSchurDel_n,
                    double *PJ, double *QJ, double *tmp_np, double *LP, double *LQ, double *Z,
                    double *u, double *w);

int cholSchurGLMvcDel(int n, int p, int r, int del_start, int del_end, double *X, double *XTilde, double *cholCap,
                      double *D1inv, double *D1invB1, double *DinvB_np, double *DinvB_nrn, double *VbetaInv,
                      double *D1inv_out, double *D1invB1_out, double *cholSchurA1_out, double *DinvB_np_out,
                      double *DinvB_nrn_out, double *cholSchurDel_n,
                      double *PB, double *QB, double *QBK, double *LP, double *LQ, double *WB, double *Z,
                      double *H, double *HK, double *A2, double *tmp_np, double *w);

void projGLMbatch(double *X, int n, int p, int b, double *V_eta, double *V_xi, double *V_beta, double *V_z,
                  double *cholpSchur, double *cholnSchur, double sigmaSqxi, double *Lbeta, double *Lz,
                  double *cholVzPlusI, double *D1invB1, double *DinvBnp, double *tmp_nb, double *tmp_pb);

void projGLMvcbatch(int n, int p, int r, int b, double *X, double *XTilde, double sigmaSqxi, double *Lbeta,
                    double *cholVz, std::string &processtype, double *V_eta, double *V_xi, double *V_beta, double *V_z,
                    double *D1inv, double *D1invB1, double *cholSchurA1_pp,
                    double *DinvB_np, double *DinvB_nrn, double *cholSchurA_nn,
                    double *tmp_nrb, double *tmp_nb, double *tmp_pb);
