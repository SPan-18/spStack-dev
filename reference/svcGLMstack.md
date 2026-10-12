# Bayesian spatially varying coefficients generalized linear model using predictive stacking

Fits Bayesian generalized linear models with spatially varying
coefficients on a collection of candidate models constructed from
candidate values of the spatial process parameters and the boundary
adjustment parameter supplied by the user, and combines inference by
stacking predictive densities. See Pan, Zhang, Bradley, and Banerjee
(2025) for more details.

## Usage

``` r
svcGLMstack(
  formula,
  data = parent.frame(),
  family,
  coords,
  cor.fn,
  process.type,
  priors,
  candidate.models,
  n.samples,
  loopd.controls,
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
  typically the environment from which `svcGLMstack` is called.

- family:

  Specifies the distribution of the response as a member of the
  exponential family. Supported options are `'poisson'`, `'binomial'`
  and `'binary'`.

- coords:

  an \\n \times 2\\ matrix of the observation coordinates in
  \\\mathbb{R}^2\\ (e.g., easting and northing).

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial dependence structure among the observations. Supported
  covariance model key words are: `'exponential'` and `'matern'`.

- process.type:

  a quoted keyword specifying the model for the spatial processes:
  `'independent'`, `'independent.shared'` or `'multivariate'`. See
  [`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md).

- priors:

  (optional) a list with each tag corresponding to a hyperparameter name
  and containing hyperprior details. Valid tags include `V.beta`,
  `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. See
  [`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md).
  If not supplied, uses defaults.

- candidate.models:

  an object of class `candidateModels` containing a list of candidate
  models for stacking, each with tags `phi`, `nu` (only if
  `cor.fn = 'matern'`), vectors of length \\r\\ if
  `process.type = 'independent'` (use
  [`list()`](https://rdrr.io/r/base/list.html) entries in
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md))
  and scalars otherwise, and `boundary`, a scalar. See
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  for details.

- n.samples:

  number of samples to be drawn from the posterior distribution.

- loopd.controls:

  a list with details on how leave-one-out predictive densities (LOO-PD)
  are to be calculated, with tags `method` (`'exact'` or `'CV'`),
  `CV.K`, `nMC` and `CV.update`, as in
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md).
  If sample size is more than 100, then the default is `'CV'` with
  `CV.K = 10`.

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

An object of class `svcGLMstack`, which is a list including the
following tags -

- `samples`:

  a list of length equal to total number of candidate models with each
  entry corresponding to a list of length 3, containing posterior
  samples of fixed effects (`beta`), the spatially varying coefficients
  (`z`), and fine scale variation `xi` for that model.

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `n.models`:

  number of candidate models that are fit.

- `model.params`:

  a list with one element per candidate model, each a named list of its
  parameters: `phi`, `nu` (`NA` for the exponential correlation
  function) and `boundary`.

- `stacking.summary`:

  a matrix with one row per candidate model, containing its parameters
  and its optimal stacking weight, for display.

- `stacking.weights`:

  a numeric vector of length equal to the number of candidate models
  storing the optimal stacking weights.

- `run.time`:

  a `proc_time` object with runtime details.

- `diagnostics`:

  a list of diagnostics, as in
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  element `numerical` (columns `min.pivot`, `min.cor` and `max.cor`, see
  [`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md))
  and element `solver`. If `verbose = TRUE`, a "Diagnostics" section is
  printed if there is an issue.

This object can be further used to recover posterior samples of the
scale parameters in the model with
[`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md),
to make predictions at new locations with
[`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md),
and to sample from the stacked posterior with
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md).

## Details

Candidate models that share the process parameters \\(\phi, \nu)\\ are
fitted together: the spatial correlation matrices and all the
pre-processing for the leave-one-out predictive densities are computed
once and shared across their boundary values, and the groups of such
candidates are the units of parallel work.

## References

Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian Inference for
Spatial-temporal Non-Gaussian Data Using Predictive Stacking." *Bayesian
Analysis*, **In Press**.
[doi:10.1214/25-BA1582](https://doi.org/10.1214/25-BA1582) .

## See also

[`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md),
[`svcLMstack()`](https://span-18.github.io/spStack-dev/reference/svcLMstack.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
data(simSpatial)
dat <- simSpatial[1:100, ]

# processes sharing (phi, nu): scalar candidates
mod.list <- candidateModels(list(phi = c(4, 8), nu = c(1.5),
                                 boundary = c(0.5, 0.75)), "cartesian")
mod1 <- svcGLMstack(y_svc_pois ~ x1 + x2 + (x1), data = dat,
                    family = "poisson",
                    coords = as.matrix(dat[, c("s1", "s2")]),
                    cor.fn = "matern", process.type = "independent.shared",
                    candidate.models = mod.list,
                    loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
                    n.samples = 500)
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | boundary | weight |
#> +---------+-----+-----+----------+--------+
#> | Model 1 |    4|  1.5|      0.50| 0.000  |
#> | Model 2 |    8|  1.5|      0.50| 0.346  |
#> | Model 3 |    4|  1.5|      0.75| 0.000  |
#> | Model 4 |    8|  1.5|      0.75| 0.654  |
#> +---------+-----+-----+----------+--------+
#> 
post_samps <- stackedSampler(mod1)
```
