#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spStack.h"

static const R_CallMethodDef CallEntries[] = {
  {"idist",        (DL_FUNC) &idist,        6},
  {"mysolveC",     (DL_FUNC) &mysolveC,     3},
  {"spGLMexact",   (DL_FUNC) &spGLMexact,   17},
  {"spLMexact",    (DL_FUNC) &spLMexact,    14},
  {"spLMexact2",   (DL_FUNC) &spLMexact2,   14},
  {"spLMexactLOO", (DL_FUNC) &spLMexactLOO, 17}
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
