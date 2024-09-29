#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);

  SEXP R_cholRankOneUpdate(SEXP L_r, SEXP n_r, SEXP v_r, SEXP alpha_r, SEXP beta_r, SEXP lower_r);

  SEXP R_cholRowDelUpdate(SEXP L_r, SEXP n_r, SEXP row_r, SEXP lower_r);

  SEXP R_cholRowBlockDelUpdate(SEXP L_r, SEXP n_r, SEXP start_r, SEXP end_r, SEXP lower_r);

  SEXP spGLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP nBinom_r,
                  SEXP coordsD_r, SEXP corfn_r, SEXP betaV_r, SEXP nu_beta_r,
                  SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP phi_r, SEXP nu_r,
                  SEXP epsilon_r, SEXP nSamples_r, SEXP verbose_r);

  SEXP spGLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP nBinom_r,
                     SEXP coordsD_r, SEXP corfn_r, SEXP betaV_r, SEXP nu_beta_r,
                     SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP phi_r, SEXP nu_r,
                     SEXP epsilon_r, SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                     SEXP CV_K_r, SEXP loopd_nMC_r, SEXP verbose_r);

  SEXP spLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                 SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                 SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                 SEXP nSamples_r, SEXP verbose_r);

  SEXP spLMexact2(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                  SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                  SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                  SEXP nSamples_r, SEXP verbose_r);

  SEXP spLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                    SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                    SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                    SEXP verbose_r);
}
