# Bayesian spatially varying coefficients linear model using predictive stacking

Fits Bayesian linear models with spatially varying coefficients for a
Gaussian response on a collection of candidate models, constructed from
candidate values of the spatial process parameters and the
noise-to-spatial variance ratios supplied by the user, and combines
inference by stacking predictive densities.

## Usage

``` r
svcLMstack(
  formula,
  data = parent.frame(),
  coords,
  cor.fn,
  process.type,
  priors = "flat",
  candidate.models,
  n.samples,
  loopd.method = "exact",
  parallel = FALSE,
  solver = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- formula:

  a symbolic description of the regression model to be fit. Variables in
  parenthesis are assigned spatially varying coefficients. See examples.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `svcLMstack` is called.

- coords:

  an \\n \times 2\\ matrix of the observation coordinates in
  \\\mathbb{R}^2\\ (e.g., easting and northing).

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial dependence structure among the observations. Supported
  covariance model key words are: `'exponential'` and `'matern'`.

- process.type:

  a quoted keyword specifying the model for the spatial processes of the
  varying coefficients: `'independent'` or `'independent.shared'`. See
  [`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md).

- priors:

  either `"flat"` (default), which assigns the prior \\p(\beta,
  \sigma^2) \propto 1/\sigma^2\\, or a list with tags `beta.norm` (a
  list containing \\\mu\_\beta\\ and \\V\_\beta\\) and/or `sigma.sq.ig`
  (a vector containing \\a\_\sigma\\ and \\b\_\sigma\\). A component not
  supplied in the list receives its flat prior.

- candidate.models:

  an object of class `candidateModels` containing a list of candidate
  models for stacking, each with tags `phi`, `nu` (only if
  `cor.fn = 'matern'`) and `noise_sp_ratio`: vectors of length \\r\\ if
  `process.type = 'independent'` (use
  [`list()`](https://rdrr.io/r/base/list.html) entries in
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)),
  otherwise scalars. See
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  for details.

- n.samples:

  number of posterior samples to be generated.

- loopd.method:

  character. Valid inputs are `'exact'` (default) and `'PSIS'`. The
  option `'exact'` finds the exact leave-one-out predictive densities in
  closed form. The option `'PSIS'` uses Pareto-smoothed importance
  sampling (Vehtari *et al.* 2024); with many latent effects its Pareto
  \\k\\ diagnostics are often high, so `'exact'` is recommended.

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

An object of class `svcLMstack`, which is a list including the following
tags -

- `samples`:

  a list of length equal to total number of candidate models with each
  entry corresponding to a list of length 4, containing posterior
  samples of fixed effects (`beta`), the noise variance (`sigmaSq`), the
  process variances (`sigmaSq.z`) and the spatially varying coefficients
  (`z`) for that model, as in
  [`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md).

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `n.models`:

  number of candidate models that are fit.

- `model.params`:

  a list with one element per candidate model, each a named list of its
  parameters: `phi`, `nu` (`NA` for the exponential correlation
  function) and `noise_sp_ratio` (vectors of length \\r\\ if
  `process.type = 'independent'`).

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
  row per candidate model (per candidate model and process, with columns
  `model` and `process`, if `process.type = 'independent'`) and columns
  `min.pivot`, `min.cor` and `max.cor`, as in
  [`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md).
  If `loopd.method = 'PSIS'`, element `pareto` contains the Pareto \\k\\
  diagnostic values of each candidate model (`k`), the threshold above
  which they are unreliable (`threshold`) and the number of values above
  it for each model (`n.high`). Element `solver` describes the
  optimization for the stacking weights: the solver used (`used`) and
  its status (`status`), the installed and requested solvers, the search
  order, the attempts with their status, and whether the fallback
  [`loo::stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  was used. If `verbose = TRUE`, a "Diagnostics" section is printed if
  there is an issue: numerical flags of candidate models with stacking
  weight above 0.05 (extreme candidates with negligible weight are
  expected in a stacking grid and are only counted), Pareto \\k\\ values
  above the threshold for any model, and solver problems.

This object can be used to make predictions at new locations with
[`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
and to sample from the stacked posterior with
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md).

## Details

Instead of assigning priors on the process parameters \\\phi\\, \\\nu\\
and the noise-to-spatial variance ratios \\\delta^2\\, we consider a set
of candidate models \\\mathcal{M} = \\M_1, \ldots, M_G\\\\ based on
candidate values of these parameters. For each \\g\\, we sample exactly
from the posterior distribution \\p(\sigma^2, \beta, z \mid y, M_g)\\
(see
[`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md))
and find the leave-one-out predictive densities \\p(y_i \mid y\_{-i},
M_g)\\. The stacking weights solve \$\$ \begin{aligned} \max\_{w_1,
\ldots, w_G}& \\ \frac{1}{n} \sum\_{i = 1}^n \log \sum\_{g = 1}^G w_g
p(y_i \mid y\_{-i}, M_g) \\ \text{subject to} & \quad w_g \geq 0,
\sum\_{g = 1}^G w_g = 1. \end{aligned} \$\$ Candidate models that share
the process parameters \\(\phi, \nu)\\ are fitted together, building and
factorizing the spatial correlation matrices only once, and the groups
of such candidates are the units of parallel work.

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

[`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md),
[`svcGLMstack()`](https://span-18.github.io/spStack-dev/reference/svcGLMstack.md),
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
[`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
data(simSpatial)
dat <- simSpatial[1:200, ]

# independent processes (r = 2): vector-valued candidates via list()
mod.list <- candidateModels(list(phi = list(c(4, 4), c(8, 4), c(8, 8)),
                                 nu = list(c(1.5, 1.5)),
                                 noise_sp_ratio = list(c(0.1, 0.1), c(0.25, 0.25))),
                            "cartesian")
mod1 <- svcLMstack(y_svc ~ x1 + x2 + (x1), data = dat,
                   coords = as.matrix(dat[, c("s1", "s2")]),
                   cor.fn = "matern", process.type = "independent",
                   candidate.models = mod.list, n.samples = 500)
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi[1] | phi[2] | nu[1] | nu[2] | noise_sp_ratio[1] | noise_sp_ratio[2] | weight |
#> +---------+--------+--------+-------+-------+-------------------+-------------------+--------+
#> | Model 1 |       4|       4|    1.5|    1.5|               0.10|               0.10| 0.000  |
#> | Model 2 |       8|       4|    1.5|    1.5|               0.10|               0.10| 0.032  |
#> | Model 3 |       8|       8|    1.5|    1.5|               0.10|               0.10| 0.206  |
#> | Model 4 |       4|       4|    1.5|    1.5|               0.25|               0.25| 0.000  |
#> | Model 5 |       8|       4|    1.5|    1.5|               0.25|               0.25| 0.762  |
#> | Model 6 |       8|       8|    1.5|    1.5|               0.25|               0.25| 0.000  |
#> +---------+--------+--------+-------+-------+-------------------+-------------------+--------+
#> 
post_samps <- stackedSampler(mod1)
n <- nrow(dat)
cor(apply(post_samps$z[n + 1:n, ], 1, median), dat$z_true)
#> [1] 0.9706991
```
