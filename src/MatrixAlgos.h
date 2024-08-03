#include <string>
#include <Rinternals.h>

void cholRankOneUpdate(int n, double *L1, double alpha, double beta,
                       double *v, double *L2, double *w);

void cholRowDelUpdate(int n, double *L, int del, double *L1, double *w);

void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1,
                 double *tmp_pp, double *tmp_np1, double*tmp_np2,
                 double *outp, double *outn);

void inversionLM2(double *X, int n, int p, double deltasq, double *VbetaInv,
                  double *Vz, double *cholVy, double *v1, double *v2,
                  double *out_p, double *out_n);

int mapIndex(int i, int j, int nRowB, int nColB, int startRowB, int startColB, int nRowA);
