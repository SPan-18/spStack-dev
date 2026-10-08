#ifndef SPSTACK_PSIS_H
#define SPSTACK_PSIS_H

int psis_tail_length(int S);

int psis_gpd_grid_length(int L);

void psis_gpdfit(double *x, int N, double *theta, double *l_theta,
                 double *k_out, double *sigma_out);

void psis_loo(double *ll, int S, int L, double *lw, int *idx, double *x_tail,
              double *theta, double *l_theta, double *elpd, double *khat);

#endif
