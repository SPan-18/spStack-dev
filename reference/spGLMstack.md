# Bayesian spatial generalized linear model using predictive stacking

Fits Bayesian spatial generalized linear model on a collection of
candidate models constructed based on some candidate values of some
model parameters specified by the user and subsequently combines
inference by stacking predictive densities. See Pan, Zhang, Bradley, and
Banerjee (2025) for more details.

## Usage

``` r
spGLMstack(
  formula,
  data = parent.frame(),
  family,
  coords,
  cor.fn,
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

  a symbolic description of the regression model to be fit. See example
  below.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `spLMstack` is called.

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
  covariance model key words are: `'exponential'` and `'matern'`. See
  below for details.

- priors:

  (optional) a list with each tag corresponding to a parameter name and
  containing prior details. Valid tags include `V.beta`, `nu.beta`,
  `nu.z` and `sigmaSq.xi`.

- candidate.models:

  an object of class `candidateModels` containing a list of candidate
  models for stacking. See
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  for details.

- n.samples:

  number of posterior samples to be generated.

- loopd.controls:

  a list with details on how leave-one-out predictive densities (LOO-PD)
  are to be calculated. Valid tags include `method`, `CV.K`, `nMC` and
  `CV.update`. The tag `method` can be either `'exact'` or `'CV'`. If
  sample size is more than 100, then the default is `'CV'` with `CV.K`
  equal to its default value 10 (Gelman *et al.* 2024). The tag `nMC`
  decides how many Monte Carlo samples will be used to evaluate the
  leave-one-out predictive densities, which must be at least 500
  (default). The tag `CV.update` is an advanced option, used only if
  `method = 'CV'`, that decides how the pre-processing of the model fit
  on each fold is obtained, and should be changed with care as the
  faster choice depends on the BLAS library that R is linked with (see
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)).
  `CV.update = 'update'` obtains it by deletion updates of the full-data
  Cholesky factors, which is the faster choice with the reference BLAS
  that R ships with. `CV.update = 'direct'` recomputes it on each fold,
  which is faster only with an optimized BLAS such as OpenBLAS, Intel
  MKL or Apple Accelerate (vecLib), and is slower otherwise. The default
  `CV.update = 'auto'` uses `'direct'` if such an optimized BLAS is
  detected from the library paths reported by R and `'update'`
  otherwise; set it explicitly if the BLAS is not detected correctly
  (for example, an optimized BLAS installed in place of `Rblas.dll` on
  Windows). Both choices give the same results up to floating-point
  rounding, so only the run time is affected.

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

An object of class `spGLMstack`, which is a list including the following
tags -

- `family`:

  the distribution of the responses as indicated in the function call

- `samples`:

  a list of length equal to total number of candidate models with each
  entry corresponding to a list of length 3, containing posterior
  samples of fixed effects (`beta`), spatial effects (`z`) and
  fine-scale variation term (`xi`) for that particular model.

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `loopd.method`:

  a list containing details of the algorithm used for calculation of
  leave-one-out predictive densities. For K-fold cross-validation, its
  tag `cv.update` records the pre-processing method (`'update'` or
  `'direct'`) that was used.

- `n.models`:

  number of candidate models that are fit.

- `model.params`:

  a list with one element per candidate model, each a named list of its
  parameters: `phi`, `nu` (`NA` for the exponential correlation
  function) and `boundary` (boundary adjustment parameter).

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
  Cholesky pivot of the \\n \times n\\ correlation matrix; values below
  1e-8 indicate a nearly singular correlation matrix), `min.cor` and
  `max.cor` (the correlations of the two farthest-apart and of the two
  closest locations; values of `min.cor` above 0.95 suggest an effective
  range far exceeding the extent of the data, values of `max.cor` below
  0.05 nearly uncorrelated locations), obtained from quantities the fit
  computes anyway. Element `solver` describes the optimization for the
  stacking weights: the solver used (`used`) and its status (`status`),
  the installed and requested solvers, the search order, the attempts
  with their status, and whether the fallback
  [`loo::stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  was used. If `verbose = TRUE`, a "Diagnostics" section is printed if
  there is an issue: numerical flags of candidate models with stacking
  weight above 0.05 (extreme candidates with negligible weight are
  expected in a stacking grid and are only counted), and solver problems
  (a requested solver not installed, an inaccurate solution, or the
  fallback).

The return object might include additional data that is useful for
subsequent prediction, model fit evaluation and other utilities.

## Details

Instead of assigning a prior on the process parameters \\\phi\\ and
\\\nu\\, the boundary adjustment parameter \\\epsilon\\, we consider a
set of candidate models based on some candidate values of these
parameters supplied by the user. Suppose the set of candidate models is
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

Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian Inference for
Spatial-temporal Non-Gaussian Data Using Predictive Stacking." *Bayesian
Analysis*, **In Press**.
[doi:10.1214/25-BA1582](https://doi.org/10.1214/25-BA1582) .

Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto Smoothed
Importance Sampling." *Journal of Machine Learning Research*,
**25**(72), 1-58. URL <https://jmlr.org/papers/v25/19-556.html>.

## See also

[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
# \donttest{
set.seed(1234)
data(simSpatial)
dat <- simSpatial[1:100, ]
cand.mod <- candidateModels(list(phi = c(3, 6, 10), nu = c(0.5, 1),
                                 boundary = c(0.5, 0.6)), "cartesian")

mod1 <- spGLMstack(y_pois ~ x1 + x2, data = dat, family = "poisson",
                   coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
                   candidate.models = cand.mod,
                   n.samples = 1000,
                   loopd.controls = list(method = "CV", CV.K = 10, nMC = 1000),
                   parallel = TRUE, verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>            | phi | nu  | boundary | weight |
#> +----------+-----+-----+----------+--------+
#> | Model 1  |    3|  0.5|       0.5| 0      |
#> | Model 2  |    6|  0.5|       0.5| 0      |
#> | Model 3  |   10|  0.5|       0.5| 0      |
#> | Model 4  |    3|  1.0|       0.5| 0      |
#> | Model 5  |    6|  1.0|       0.5| 0      |
#> | Model 6  |   10|  1.0|       0.5| 0      |
#> | Model 7  |    3|  0.5|       0.6| 0      |
#> | Model 8  |    6|  0.5|       0.6| 0      |
#> | Model 9  |   10|  0.5|       0.6| 0      |
#> | Model 10 |    3|  1.0|       0.6| 0      |
#> | Model 11 |    6|  1.0|       0.6| 0      |
#> | Model 12 |   10|  1.0|       0.6| 1      |
#> +----------+-----+-----+----------+--------+
#> 

post_samps <- stackedSampler(mod1)
post_beta <- post_samps$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>                   2.5%        50%      97.5%
#> (Intercept)  0.5963552  1.7859907  2.7538797
#> x1          -0.7582803 -0.5239114 -0.3066513
#> x2           0.1323709  0.4129337  0.7259990

# compare the posterior medians of the spatial effects with the truth
cor(apply(post_samps$z, 1, median), dat$z_true)
#> [1] 0.9179889
# }
```
