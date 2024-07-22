#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);

  SEXP mysolveC(SEXP A_r, SEXP b_r, SEXP n_r);

  SEXP spLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                 SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                 SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                 SEXP nSamples_r, SEXP verbose_r);

  SEXP spLMexact2(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                  SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                  SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                  SEXP nSamples_r, SEXP verbose_r);
}
