#include <string>
#include <Rinternals.h>

// SEXP mysolve(SEXP A_r, SEXP b_r);
void mysolve(double *A, double *b, int n);
