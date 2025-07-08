#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);

  SEXP recoverScale_stvcGLM(SEXP n_r, SEXP p_r, SEXP r_r, SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                            SEXP betaMu_r, SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP iwScale_r, SEXP processType_r,
                            SEXP phi_s_r, SEXP phi_t_r, SEXP nSamples_r, SEXP betaSamps_r, SEXP zSamps_r);

  SEXP recoverScale_spGLM(SEXP n_r, SEXP p_r, SEXP coords_r, SEXP corfn_r,
                          SEXP betaMu_r, SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r,
                          SEXP phi_r, SEXP nu_r, SEXP nSamples_r, SEXP betaSamps_r, SEXP zSamps_r);

  SEXP R_cholRankOneUpdate(SEXP L_r, SEXP n_r, SEXP v_r, SEXP alpha_r, SEXP beta_r, SEXP lower_r);

  SEXP R_cholRowDelUpdate(SEXP L_r, SEXP n_r, SEXP row_r, SEXP lower_r);

  SEXP R_cholRowBlockDelUpdate(SEXP L_r, SEXP n_r, SEXP start_r, SEXP end_r, SEXP lower_r);

  SEXP predict_spGLM(SEXP n_r, SEXP n_pred_r, SEXP p_r, SEXP family_r, SEXP nBinom_new_r,
                     SEXP X_new_r, SEXP sp_coords_r, SEXP sp_coords_new_r,
                     SEXP corfn_r, SEXP phi_r, SEXP nu_r, SEXP nSamples_r,
                     SEXP beta_samps_r, SEXP z_samps_r, SEXP sigmaSq_z_samps_r, SEXP joint_r);

  SEXP predict_stvcGLM(SEXP n_r, SEXP n_pred_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_new_r, SEXP X_new_r, SEXP XTilde_new_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP sp_coords_new_r, SEXP time_coords_new_r,
                       SEXP processType_r, SEXP corfn_r, SEXP phi_s_r, SEXP phi_t_r, SEXP nSamples_r,
                       SEXP beta_samps_r, SEXP z_samps_r, SEXP z_scale_samps_r, SEXP joint_r);

  SEXP predict_spLM(SEXP n_r, SEXP n_pred_r, SEXP p_r,
                    SEXP X_new_r, SEXP sp_coords_r, SEXP sp_coords_new_r,
                    SEXP corfn_r, SEXP phi_r, SEXP nu_r, SEXP deltasq_r,
                    SEXP beta_samps_r, SEXP z_samps_r, SEXP sigmaSq_z_samps_r,
                    SEXP nSamples_r, SEXP joint_r);

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

  SEXP stvcGLMexact(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                    SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                    SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                    SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                    SEXP nSamples_r, SEXP verbose_r);

  SEXP stvcGLMexactLOO(SEXP Y_r, SEXP X_r, SEXP X_tilde_r, SEXP n_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP nBinom_r,
                       SEXP sp_coords_r, SEXP time_coords_r, SEXP corfn_r,
                       SEXP betaV_r, SEXP nu_beta_r, SEXP nu_z_r, SEXP sigmaSq_xi_r, SEXP iwScale_r,
                       SEXP processType_r, SEXP phi_s_r, SEXP phi_t_r, SEXP epsilon_r,
                       SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                       SEXP CV_K_r, SEXP loopd_nMC_r,  SEXP verbose_r);
}
