#include <string>
#include <Rinternals.h>

void copyMatrixDelRow(double *M1, int nRowM1, int nColM1, double *M2, int exclude_index);

void copyMatrixColDelRowBlock(double *M1, int nRowM1, int nColM1, double *M2,
                              int include_start, int include_end, int exclude_start, int exclude_end);

void copyMatrixDelRowBlock(double *M1, int nRowM1, int nColM1, double *M2, int exclude_start, int exclude_end);

void copyMatrixDelRowCol(double *M1, int nRowM1, int nColM1, double *M2, int del_indexRow, int del_indexCol);

void copyMatrixDelRowColBlock(double *M1, int nRowM1, int nColM1, double *M2,
                              int delRow_start, int delRow_end, int delCol_start, int delCol_end);

void copyMatrixRowBlock(double *M1, int nRowM1, int nColM1, double *M2, int copy_start, int copy_end);

void copyMatrixDelRowBlock_vc(double *M1, int nRowM1, int nColM1, double *M2, int exclude_start, int exclude_end, int rep);

void copyMatrixRowColBlock(double *M1, int nRowM1, int nColM1, double *M2,
                           int copyCol_start, int copyCol_end, int copyRow_start, int copyRow_end);

void copyMatrixDelRowColBlock_vc(double *M1, int nRowM1, int nColM1, double *M2, int delRow_start, int delRow_end,
                                 int delCol_start, int delCol_end, int rep);

void copyMatrixRowToVec(double *M, int nRowM, int nColM, double *vec, int copy_index);

void copySubmat(double *A, int nRowA, int nColA, double *B, int nRowB, int nColB,
                int startRowA, int startColA, int startRowB, int startColB,
                int nRowCopy, int nColCopy);

void copyVecBlock(double *v1, double *v2, int n, int copy_start, int copy_end);

void copyVecExcludingBlock(double *v1, double *v2, int n, int exclude_start, int exclude_end);

void copyVecExcludingOne(double *v1, double *v2, int n, int exclude_index);

int findMax(int *a, int n);

double findMax(double *a, int n);

int findMin(int *a, int n);

double inverse_logit(double x);

double logit(double x);

double rlogGamma(double shape);

double rlogitBeta(double a, double b);

double logMeanExp(double *a, int n);

double logSumExp(double *a, int n);

void mkCVpartition(int n, int K, int *start_vec, int *end_vec, int *size_vec);

void mkLT(double *A, int n);

void printMtrx(double *m, int nRow, int nCol);

void printVec(double *m, int n);

void printVec(int *m, int n);

void spCorFull2(int n, int p, double *coords_sp, double *theta, std::string &corfn, double *C);

void sptCorFull(int n, int p, double *coords_sp, double *coords_tm, double *theta, std::string &corfn, double *C);

void spCorCross(int n, int n_prime, int p, double *coords_sp, double *coords_sp_prime, double *theta, std::string &corfn, double *C);

void sptCorCross(int n, int n_prime, int p, double *coords_sp, double *coords_tm,
                 double *coords_sp_prime, double *coords_tm_prime, double *theta, std::string &corfn, double *C);

double gneiting_spt_decay(double dist_s, double dist_t, double phi_s, double phi_t);

void zeros(double *x, int length);

void zeros(int *x, int length);

void lmulv_XTilde_VC(const char *trans, int n, int r, double *XTilde, double *v, double *res);

void lmulm_XTilde_VC(const char *trans, int n, int r, int k, double *XTilde, double *A, double *res);

void rmul_Vz_XTildeT(int n, int r, double *XTilde, double *Vz, double *res, std::string &processtype);

void rWishartBartlett(int r, double nu, double *A);

int invWishartFromBartlett(int r, double *A, double *cholinvIWscale, double *Sigma, double *tmp_rr);

int rInvWishart(int r, double nu, double *cholinvIWscale, double *Sigma, double *tmp_rr);
void corOffDiagRange(double *A, int n, double *minCor, double *maxCor);

double minRelPivot(double *L, int n, double *d, double dconst);

SEXP appendDiagnostics(SEXP list_r, double minPivot, double minCor, double maxCor);

SEXP appendDiagnosticsRows(SEXP list_r, double *minPivot, double *minCor, double *maxCor, int m);

int pendingInterrupt();
