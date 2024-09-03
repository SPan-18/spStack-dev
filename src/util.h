#include <string>
#include <Rinternals.h>

void copyVecExcludingOne(double *v1, double *v2, int n, int exclude_index);

void copyMatrixDelRow(double *M1, int nRowM1, int nColM1, double *M2, int exclude_index);

void copyMatrixDelRowCol(double *M1, int nRowM1, int nColM1, double *M2, int del_indexRow, int del_indexCol);

void copyMatrixRowToVec(double *M, int nRowM, int nColM, double *vec, int copy_index);

void copyMatrixSEXP(double *matrixC, int dim1, int dim2, double *pointerSEXP);

void copySubmat(double *A, int nRowA, int nColA, double *B, int nRowB, int nColB,
                int startRowA, int startColA, int startRowB, int startColB,
                int nRowCopy, int nColCopy);

void copyVectorSEXP(double *vectorC, int dim, double *pointerSEXP);

double logit(double x);

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
