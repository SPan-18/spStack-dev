#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "spStack.h"

static const R_CallMethodDef CallEntries[] = {
  {"idist",                   (DL_FUNC) &idist,                   6},
  {"recoverScale_stvcGLM",    (DL_FUNC) &recoverScale_stvcGLM,    17},
  {"recoverScale_spGLM",      (DL_FUNC) &recoverScale_spGLM,      13},
  {"R_cholRankOneUpdate",     (DL_FUNC) &R_cholRankOneUpdate,     6},
  {"R_cholRowDelUpdate",      (DL_FUNC) &R_cholRowDelUpdate,      4},
  {"R_cholRowBlockDelUpdate", (DL_FUNC) &R_cholRowBlockDelUpdate, 5},
  {"R_psis",                  (DL_FUNC) &R_psis,                  2},
  {"predict_spGLM",           (DL_FUNC) &predict_spGLM,           16},
  {"predict_stvcGLM",         (DL_FUNC) &predict_stvcGLM,         21},
  {"predict_stvcLM",          (DL_FUNC) &predict_stvcLM,          20},
  {"predict_spLM",            (DL_FUNC) &predict_spLM,            15},
  {"spGLMexact",              (DL_FUNC) &spGLMexact,              17},
  {"spGLMexactLOO",           (DL_FUNC) &spGLMexactLOO,           22},
  {"spGLMexactLOOgrid",       (DL_FUNC) &spGLMexactLOOgrid,       21},
  {"spLMexact",               (DL_FUNC) &spLMexact,               14},
  {"spLMexactLOO",            (DL_FUNC) &spLMexactLOO,            16},
  {"spLMexactLOOgrid",        (DL_FUNC) &spLMexactLOOgrid,        15},
  {"stvcGLMexact",            (DL_FUNC) &stvcGLMexact,            22},
  {"stvcGLMexactLOO",         (DL_FUNC) &stvcGLMexactLOO,         27},
  {"stvcGLMexactLOOgrid",     (DL_FUNC) &stvcGLMexactLOOgrid,     27},
  {"stvcLMexact",             (DL_FUNC) &stvcLMexact,             20},
  {"stvcLMexactGrid",         (DL_FUNC) &stvcLMexactGrid,         19},
  {NULL, NULL, 0}
};

// Entry point called by R when the shared library is loaded
extern "C" void attribute_visible R_init_spStack(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
