#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);

  SEXP mysolveC(SEXP A_r, SEXP b_r, SEXP n_r);

}
