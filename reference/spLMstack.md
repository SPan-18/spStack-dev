# Bayesian spatial linear model using predictive stacking

Fits Bayesian spatial linear model on a collection of candidate models
constructed based on some candidate values of some model parameters
specified by the user and subsequently combines inference by stacking
predictive densities. See Zhang, Tang and Banerjee (2025) for more
details.

## Usage

``` r
spLMstack(
  formula,
  data = parent.frame(),
  coords,
  cor.fn,
  priors,
  params.list,
  n.samples,
  loopd.method,
  parallel = FALSE,
  solver = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- formula:

  a symbolic description of the regression model to be fit. See example
  below.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `spLMstack` is called.

- coords:

  an \\n \times 2\\ matrix of the observation coordinates in
  \\\mathbb{R}^2\\ (e.g., easting and northing).

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial dependence structure among the observations. Supported
  covariance model key words are: `'exponential'` and `'matern'`. See
  below for details.

- priors:

  a list with each tag corresponding to a parameter name and containing
  prior details. If not supplied, uses defaults.

- params.list:

  a list containing candidate values of spatial process parameters for
  the `cor.fn` used, and, noise-to-spatial variance ratio.

- n.samples:

  number of posterior samples to be generated.

- loopd.method:

  character. Valid inputs are `'exact'` and `'PSIS'`. The option
  `'exact'` corresponds to exact leave-one-out predictive densities. The
  option `'PSIS'` is faster, as it finds approximate leave-one-out
  predictive densities using Pareto-smoothed importance sampling (Gelman
  *et al.* 2024).

- parallel:

  logical. If `parallel=FALSE`, the parallelization plan, if set up by
  the user, is ignored. If `parallel=TRUE`, the function inherits the
  parallelization plan that is set by the user via the function
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  only. Depending on the parallel backend available, users may choose
  their own plan. More details are available at
  <https://cran.R-project.org/package=future>.

- solver:

  (optional) Specifies the name of the solver that will be used to
  obtain optimal stacking weights for each candidate model. Default
  order is `c("CLARABEL", "ECOS", "SCS")`. Users can use other solvers
  supported by the
  [CVXR-package](https://www.cvxgrp.org/CVXR/reference/CVXR-package.html)
  package.

- verbose:

  logical. If `TRUE`, prints model-specific optimal stacking weights.

- ...:

  currently no additional argument.

## Value

An object of class `spLMstack`, which is a list including the following
tags -

- `samples`:

  a list of length equal to total number of candidate models with each
  entry corresponding to a list of length 3, containing posterior
  samples of fixed effects (`beta`), variance parameter (`sigmaSq`),
  spatial effects (`z`) for that model.

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `n.models`:

  number of candidate models that are fit.

- `candidate.models`:

  a matrix with `n_model` rows with each row containing details of the
  model parameters and its optimal weight.

- `stacking.weights`:

  a numeric vector of length equal to the number of candidate models
  storing the optimal stacking weights.

- `run.time`:

  a `proc_time` object with runtime details.

- `solver.status`:

  solver status as returned by the optimization routine.

The return object might include additional data that is useful for
subsequent prediction, model fit evaluation and other utilities.

## Details

Instead of assigning a prior on the process parameters \\\phi\\ and
\\\nu\\, noise-to-spatial variance ratio \\\delta^2\\, we consider a set
of candidate models based on some candidate values of these parameters
supplied by the user. Suppose the set of candidate models is
\\\mathcal{M} = \\M_1, \ldots, M_G\\\\. Then for each \\g = 1, \ldots,
G\\, we sample from the posterior distribution \\p(\sigma^2, \beta, z
\mid y, M_g)\\ under the model \\M_g\\ and find leave-one-out predictive
densities \\p(y_i \mid y\_{-i}, M_g)\\. Then we solve the optimization
problem \$\$ \begin{aligned} \max\_{w_1, \ldots, w_G}& \\ \frac{1}{n}
\sum\_{i = 1}^n \log \sum\_{g = 1}^G w_g p(y_i \mid y\_{-i}, M_g) \\
\text{subject to} & \quad w_g \geq 0, \sum\_{g = 1}^G w_g = 1
\end{aligned} \$\$ to find the optimal stacking weights \\\hat{w}\_1,
\ldots, \hat{w}\_G\\.

## References

Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto Smoothed
Importance Sampling." *Journal of Machine Learning Research*,
**25**(72), 1-58. URL <https://jmlr.org/papers/v25/19-556.html>.

Zhang L, Tang W, Banerjee S (2025). "Bayesian Geostatistics Using
Predictive Stacking." *Journal of the American Statistical Association*,
**In press**.
[doi:10.1080/01621459.2025.2566449](https://doi.org/10.1080/01621459.2025.2566449)
.

## See also

[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
# load data and work with first 100 rows
data(simGaussian)
dat <- simGaussian[1:100, ]

# setup prior list
muBeta <- c(0, 0)
VBeta <- cbind(c(1.0, 0.0), c(0.0, 1.0))
sigmaSqIGa <- 2
sigmaSqIGb <- 2
prior_list <- list(beta.norm = list(muBeta, VBeta),
                   sigma.sq.ig = c(sigmaSqIGa, sigmaSqIGb))

mod1 <- spLMstack(y ~ x1, data = dat,
                  coords = as.matrix(dat[, c("s1", "s2")]),
                  cor.fn = "matern",
                  priors = prior_list,
                  params.list = list(phi = c(1.5, 3),
                                     nu = c(0.5, 1),
                                     noise_sp_ratio = c(1)),
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.8.1 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.05s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -29.233
#> ℹ Compile time: 0.05s
#> ℹ Solver time: 0.005s
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | noise_sp_ratio | weight |
#> +---------+-----+-----+----------------+--------+
#> | Model 1 |  1.5|  0.5|               1| 0      |
#> | Model 2 |  3.0|  0.5|               1| 0      |
#> | Model 3 |  1.5|  1.0|               1| 0      |
#> | Model 4 |  3.0|  1.0|               1| 1      |
#> +---------+-----+-----+----------------+--------+
#> 

post_samps <- stackedSampler(mod1)
post_beta <- post_samps$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>                 2.5%      50%    97.5%
#> (Intercept) 0.878944 1.739341 2.612417
#> x1          4.745753 4.925319 5.076682

post_z <- post_samps$z
post_z_summ <- t(apply(post_z, 1,
                       function(x) quantile(x, c(0.025, 0.5, 0.975))))

z_combn <- data.frame(z = dat$z_true,
                      zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2],
                      zU = post_z_summ[, 3])

library(ggplot2)
plot1 <- ggplot(data = z_combn, aes(x = z)) +
  geom_point(aes(y = zM), size = 0.25,
             color = "darkblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = zL, ymax = zU),
                width = 0.05, alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "solid") +
  xlab("True z") + ylab("Stacked posterior of z") +
  theme_bw() +
  theme(panel.background = element_blank(),
        aspect.ratio = 1)
```
