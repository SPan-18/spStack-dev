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

void copyMatrixRowColBlock(double *M1, int nRowM1, int nColM1, double *M2,
                           int copyCol_start, int copyCol_end, int copyRow_start, int copyRow_end);

void copyMatrixRowToVec(double *M, int nRowM, int nColM, double *vec, int copy_index);

void copyMatrixSEXP(double *matrixC, int dim1, int dim2, double *pointerSEXP);

void copySubmat(double *A, int nRowA, int nColA, double *B, int nRowB, int nColB,
                int startRowA, int startColA, int startRowB, int startColB,
                int nRowCopy, int nColCopy);

void copyVecBlock(double *v1, double *v2, int n, int copy_start, int copy_end);

void copyVecExcludingBlock(double *v1, double *v2, int n, int exclude_start, int exclude_end);

void copyVecExcludingOne(double *v1, double *v2, int n, int exclude_index);

void copyVectorSEXP(double *vectorC, int dim, double *pointerSEXP);

int findMax(int *a, int n);

int findMin(int *a, int n);

double inverse_logit(double x);

double logit(double x);

double logMeanExp(double *a, int n);

void mkCVpartition(int n, int K, int *start_vec, int *end_vec, int *size_vec);

void mkLT(double *A, int n);

void mysolveLT(double *A, double *b, int n);

void mysolveUT(double *A, double *b, int n);

void printMtrx(double *m, int nRow, int nCol);

void printVec(double *m, int n);

void printVec(int *m, int n);

void spCorLT(double *D, int n, double *theta, std::string &corfn, double *C);

void spCorFull(double *D, int n, double *theta, std::string &corfn, double *C);

void zeros(double *x, int length);

void zeros(int *x, int length);
