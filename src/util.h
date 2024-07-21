#include <string>
#include <Rinternals.h>

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
