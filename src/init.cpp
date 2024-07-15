#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spStack.h"

static const R_CallMethodDef CallEntries[] = {
  {"idist",      (DL_FUNC) &idist,       6},
  {"mysolveC",   (DL_FUNC) &mysolveC,    3},
  {"spLM_fixed", (DL_FUNC) &spLM_fixed, 14}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
  R_init_sp(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
  }
