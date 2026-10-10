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
  priors = "flat",
  candidate.models,
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

  either `"flat"` (default), which assigns the prior \\p(\beta,
  \sigma^2) \propto 1/\sigma^2\\, or a list with tags `beta.norm` (a
  list containing \\\mu\_\beta\\ and \\V\_\beta\\) and/or `sigma.sq.ig`
  (a vector containing \\a\_\sigma\\ and \\b\_\sigma\\). A component not
  supplied in the list receives its flat prior, \\p(\beta) \propto 1\\
  or \\p(\sigma^2) \propto 1/\sigma^2\\.

- candidate.models:

  an object of class `candidateModels` containing a list of candidate
  models for stacking. See
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  for details.

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
  entry corresponding to a list of length 4, containing posterior
  samples of fixed effects (`beta`), measurement error variance
  (`sigmaSq`), spatial variance (`sigmaSq.z`), and spatial effects (`z`)
  for that model.

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `n.models`:

  number of candidate models that are fit.

- `model.params`:

  a list with one element per candidate model, each a named list of its
  parameters: `phi`, `nu` (`NA` for the exponential correlation
  function) and `noise_sp_ratio` (noise-to-spatial variance ratio).

- `stacking.summary`:

  a matrix with one row per candidate model, containing its parameters
  and its optimal stacking weight, for display.

- `stacking.weights`:

  a numeric vector of length equal to the number of candidate models
  storing the optimal stacking weights.

- `run.time`:

  a `proc_time` object with runtime details.

- `diagnostics`:

  a list of diagnostics. Element `numerical` is a data frame with one
  row per candidate model and columns `min.pivot` (the smallest relative
  Cholesky pivot of the \\n \times n\\ factorizations; values below 1e-8
  indicate a nearly singular covariance matrix), `min.cor` and `max.cor`
  (the correlations of the two farthest-apart and of the two closest
  locations; values of `min.cor` above 0.95 suggest an effective range
  far exceeding the extent of the data, values of `max.cor` below 0.05
  nearly uncorrelated locations), obtained from quantities the fits
  compute anyway. If `loopd.method = 'PSIS'`, element `pareto` is a list
  with the Pareto \\k\\ diagnostic values of each candidate model (`k`),
  the threshold above which they are unreliable (`threshold`) and the
  number of values above it for each model (`n.high`). Element `solver`
  describes the optimization for the stacking weights: the solver used
  (`used`) and its status (`status`), the installed and requested
  solvers, the search order, the attempts with their status, and whether
  the fallback
  [`loo::stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  was used. If `verbose = TRUE`, a "Diagnostics" section is printed if
  there is an issue: numerical flags of candidate models with stacking
  weight above 0.05 (extreme candidates with negligible weight are
  expected in a stacking grid and are only counted), Pareto \\k\\ values
  above the threshold for any model, and solver problems (a requested
  solver not installed, an inaccurate solution, or the fallback).

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
data(simSpatial)
dat <- simSpatial[1:100, ]

cand.mod <- candidateModels(list(phi = c(3, 6),
                                 nu = c(0.5, 1),
                                 noise_sp_ratio = c(0.5, 1)),
                            "cartesian")

mod1 <- spLMstack(y_gauss ~ x1 + x2, data = dat,
                  coords = as.matrix(dat[, c("s1", "s2")]),
                  cor.fn = "matern",
                  candidate.models = cand.mod,
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | noise_sp_ratio | weight |
#> +---------+-----+-----+----------------+--------+
#> | Model 1 |    3|  0.5|             0.5| 0      |
#> | Model 2 |    6|  0.5|             0.5| 0      |
#> | Model 3 |    3|  1.0|             0.5| 0      |
#> | Model 4 |    6|  1.0|             0.5| 1      |
#> | Model 5 |    3|  0.5|             1.0| 0      |
#> | Model 6 |    6|  0.5|             1.0| 0      |
#> | Model 7 |    3|  1.0|             1.0| 0      |
#> | Model 8 |    6|  1.0|             1.0| 0      |
#> +---------+-----+-----+----------------+--------+
#> 

post_samps <- stackedSampler(mod1)
post_beta <- post_samps$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>                  2.5%        50%      97.5%
#> (Intercept)  1.528014  2.2337693  2.9973924
#> x1           4.729069  4.8543561  4.9913745
#> x2          -1.130784 -0.9844785 -0.8387339

# compare the posterior medians of the spatial effects with the truth
cor(apply(post_samps$z, 1, median), dat$z_true)
#> [1] 0.9574472
```
