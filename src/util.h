#include <string>
#include <Rinternals.h>

// SEXP mysolve(SEXP A_r, SEXP b_r);
void mysolve(double *A, double *b, int n);

void printMtrx(double *m, int nRow, int nCol);

void printVec(double *m, int n);

void printVec(int *m, int n);
