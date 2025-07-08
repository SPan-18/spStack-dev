#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spStack.h"

static const R_CallMethodDef CallEntries[] = {
  {"idist",                   (DL_FUNC) &idist,                   6},
  {"recoverScale_stvcGLM",    (DL_FUNC) &recoverScale_stvcGLM,    17},
  {"recoverScale_spGLM",      (DL_FUNC) &recoverScale_spGLM,      13},
  {"R_cholRankOneUpdate",     (DL_FUNC) &R_cholRankOneUpdate,     6},
  {"R_cholRowDelUpdate",      (DL_FUNC) &R_cholRowDelUpdate,      4},
  {"R_cholRowBlockDelUpdate", (DL_FUNC) &R_cholRowBlockDelUpdate, 5},
  {"predict_spGLM",           (DL_FUNC) &predict_spGLM,           16},
  {"predict_stvcGLM",         (DL_FUNC) &predict_stvcGLM,         21},
  {"predict_spLM",            (DL_FUNC) &predict_spLM,            15},
  {"spGLMexact",              (DL_FUNC) &spGLMexact,              17},
  {"spGLMexactLOO",           (DL_FUNC) &spGLMexactLOO,           21},
  {"spLMexact",               (DL_FUNC) &spLMexact,               14},
  {"spLMexact2",              (DL_FUNC) &spLMexact2,              14},
  {"spLMexactLOO",            (DL_FUNC) &spLMexactLOO,            16},
  {"stvcGLMexact",            (DL_FUNC) &stvcGLMexact,            22},
  {"stvcGLMexactLOO",         (DL_FUNC) &stvcGLMexactLOO,         26}
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
