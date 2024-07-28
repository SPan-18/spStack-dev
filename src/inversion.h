#include <string>
#include <Rinternals.h>

void inversionLM(double *X, int n, int p, double deltasq, double *VbetaInv,
                 double *Vz, double *cholVy, double *v1, double *v2,
                 double *tmp_n1, double *tmp_n2, double *tmp_p1,
                 double *tmp_pp, double *tmp_np1, double*tmp_np2,
                 double *outp, double *outn);

void inversionLM2(double *X, int n, int p, double deltasq, double *VbetaInv,
                  double *Vz, double *cholVy, double *v1, double *v2,
                  double *out_p, double *out_n);

void zeros2(double *x, int length);
