#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include "MatrixAlgos.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP predict_stvcGLMexact(SEXP n_r, SEXP n_pred_r, SEXP p_r, SEXP r_r, SEXP family_r, SEXP X_new_r, SEXP XTilde_new_r,
                           SEXP sp_coords_r, SEXP time_coords_r, SEXP sp_coords_new_r, SEXP time_coords_new_r,
                           SEXP processType_r, SEXP corfn_r, SEXP phi_s_r, SEXP phi_t_r, SEXP nSamples_r,
                           SEXP beta_samps_r, SEXP z_samps_r, SEXP sigmaSqbeta_samps_r, SEXP z_scale_samps_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, k, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
     Set-up
     *****************************************/
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int n_pred = INTEGER(n_pred_r)[0];
    int n_predn_pred = n_pred * n_pred;
    int nn_pred = n * n_pred;
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int r = INTEGER(r_r)[0];
    int rr = r * r;
    int nr = n * r;
    int nnr = nn * r;
    int n_predn_predr = n_predn_pred * r;
    int nn_predr = nn_pred * r;

    double *zSamps = REAL(z_samps_r);
    double *betaSamps = REAL(beta_samps_r);

    double *coords_sp = REAL(sp_coords_r);
    double *coords_sp_new = REAL(sp_coords_new_r);
    double *coords_tm = REAL(time_coords_r);
    double *coords_tm_new = REAL(time_coords_new_r);

    std::string corfn = CHAR(STRING_ELT(corfn_r, 0));

    // create spatial-temporal covariance matrices
    std::string processType = CHAR(STRING_ELT(processType_r, 0));
    double *phi_s_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_s_vec, r);
    double *phi_t_vec = (double *) R_alloc(r, sizeof(double)); zeros(phi_t_vec, r);
    double *thetaspt = (double *) R_alloc(2, sizeof(double));
    double *Vz = NULL;
    double *Vz_new = NULL;
    double *Cz = NULL;

    if(corfn == "gneiting-decay"){

        if(processType == "independent.shared" || processType == "multivariate"){

            phi_s_vec[0] = REAL(phi_s_r)[0];
            phi_t_vec[0] = REAL(phi_t_r)[0];
            thetaspt[0] = phi_s_vec[0];
            thetaspt[1] = phi_t_vec[0];

            Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);
            Vz_new = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(Vz_new, n_predn_pred);
            Cz = (double *) R_alloc(nn_pred, sizeof(double)); zeros(Cz, nn_pred);

            sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, Vz);
            sptCorFull(n_pred, 2, coords_sp_new, coords_tm_new, thetaspt, corfn, Vz_new);
            sptCorCross(n, n_pred, 2, coords_sp, coords_tm, coords_sp_new, coords_tm_new, thetaspt, corfn, Cz);

        }else if(processType == "independent"){

            F77_NAME(dcopy)(&r, REAL(phi_s_r), &incOne, phi_s_vec, &incOne);
            F77_NAME(dcopy)(&r, REAL(phi_t_r), &incOne, phi_t_vec, &incOne);

            Vz = (double *) R_alloc(nnr, sizeof(double)); zeros(Vz, nnr);
            Vz_new = (double *) R_alloc(n_predn_predr, sizeof(double)); zeros(Vz_new, n_predn_predr);
            Cz = (double *) R_alloc(nn_predr, sizeof(double)); zeros(Cz, nn_predr);

            // find r-many correlation/cross-correlation matrices, stacked into a rn^2-dim vector
            for(k = 0; k < r; k++){
                thetaspt[0] = phi_s_vec[k];
                thetaspt[1] = phi_t_vec[k];
                sptCorFull(n, 2, coords_sp, coords_tm, thetaspt, corfn, &Vz[nn * k]);
                sptCorFull(n_pred, 2, coords_sp_new, coords_tm_new, thetaspt, corfn, &Vz_new[n_predn_pred * k]);
                sptCorCross(n, n_pred, 2, coords_sp, coords_tm, coords_sp_new, coords_tm_new, thetaspt, corfn, &Cz[nn_pred * k]);
            }
        }
    }

    // sampling set-up
    int nSamples = INTEGER(nSamples_r)[0];

    /*****************************************
     Set-up preprocessing matrices etc.
     *****************************************/

    double *cholVz = NULL;               // define NULL pointer for chol(Vz)
    double *cholVz_new = NULL;           // define NULL pointer for chol(Vz_new)
    double *z_pred_cov = NULL;           // define NULL pointer for z_pred_cov
    double *z_pred_mu = NULL;            // define NULL pointer for z_pred_mu

    // Find Cholesky of Vz
    if(processType == "independent.shared"){

        cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);                              // nxn matrix chol(Vz)
        cholVz_new = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(cholVz_new, n_predn_pred);  // n_predxn_pred matrix chol(Vz_new)

        F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dcopy)(&n_predn_pred, Vz_new, &incOne, cholVz_new, &n_predn_pred);

        F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
        F77_NAME(dpotrf)(lower, &n_pred, cholVz_new, &n_pred, &info FCONE); if(info != 0){perror("c++ error: Vz_new dpotrf failed\n");}
        mkLT(cholVz, n);
        mkLT(cholVz_new, n_pred);

    }else if(processType == "independent"){

        cholVz = (double *) R_alloc(nnr, sizeof(double)); zeros(cholVz, nnr);                             // r nxn matrices chol(Vz)
        cholVz_new = (double *) R_alloc(n_predn_predr, sizeof(double)); zeros(cholVz_new, n_predn_predr); // r n_predxn_pred matrices chol(Vz_new)

        F77_NAME(dcopy)(&nnr, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dcopy)(&n_predn_predr, Vz_new, &incOne, cholVz_new, &incOne);

        for(k = 0; k < r; k++){
            F77_NAME(dpotrf)(lower, &n, &cholVz[nn * k], &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
            F77_NAME(dpotrf)(lower, &n_pred, &cholVz_new[n_predn_pred * k], &n_pred, &info FCONE); if(info != 0){perror("c++ error: Vz_new dpotrf failed\n");}
            mkLT(&cholVz[nn * k], n);
            mkLT(&cholVz_new[n_predn_pred * k], n_pred);
        }

    }else if(processType == "multivariate"){

        cholVz = (double *) R_alloc(nn, sizeof(double)); zeros(cholVz, nn);                              // nxn matrix chol(Vz)
        cholVz_new = (double *) R_alloc(n_predn_pred, sizeof(double)); zeros(cholVz_new, n_predn_pred);  // n_predxn_pred matrix chol(Vz_new)

        F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
        F77_NAME(dcopy)(&n_predn_pred, Vz_new, &incOne, cholVz_new, &n_predn_pred);

        F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
        F77_NAME(dpotrf)(lower, &n_pred, cholVz_new, &n_pred, &info FCONE); if(info != 0){perror("c++ error: Vz_new dpotrf failed\n");}
        mkLT(cholVz, n);
        mkLT(cholVz_new, n_pred);

    }

    return R_NilValue;

  }

}