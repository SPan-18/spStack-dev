#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP mysolveC(SEXP A_r, SEXP b_r, SEXP n_r);

}
