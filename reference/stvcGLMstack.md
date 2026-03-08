# Bayesian spatially-temporally varying coefficients generalized linear model using predictive stacking

Fits Bayesian spatial-temporal generalized linear model with
spatially-temporally varying coefficients on a collection of candidate
models constructed based on some candidate values of some model
parameters specified by the user and subsequently combines inference by
stacking predictive densities. See Pan, Zhang, Bradley, and Banerjee
(2025) for more details.

## Usage

``` r
stvcGLMstack(
  formula,
  data = parent.frame(),
  family,
  sp_coords,
  time_coords,
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
  parenthesis are assigned spatially-temporally varying coefficients.
  See examples.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `stvcGLMstack` is called.

- family:

  Specifies the distribution of the response as a member of the
  exponential family. Supported options are `'poisson'`, `'binomial'`
  and `'binary'`.

- sp_coords:

  an \\n \times 2\\ matrix of the observation spatial coordinates in
  \\\mathbb{R}^2\\ (e.g., easting and northing).

- time_coords:

  an \\n \times 1\\ matrix of the observation temporal coordinates in
  \\\mathcal{T} \subseteq \[0, \infty)\\.

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial-temporal dependence structure among the observations.
  Supported covariance model key words are: `'gneiting-decay'` (Gneiting
  and Guttorp 2010). See below for details.

- process.type:

  a quoted keyword specifying the model for the spatial-temporal
  process. Supported keywords are `'independent'` which indicates
  independent processes for each varying coefficients characterized by
  different process parameters, `independent.shared` implies independent
  processes for the varying coefficients that shares common process
  parameters, and `multivariate` implies correlated processes for the
  varying coefficients modeled by a multivariate Gaussian process with
  an inverse-Wishart prior on the correlation matrix. The input for
  `sptParams` and `priors` must be given accordingly.

- priors:

  (optional) a list with each tag corresponding to a hyperparameter name
  and containing hyperprior details. Valid tags include `V.beta`,
  `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
  and `nu.z` must be at least 2.1. If not supplied, uses defaults.

- candidate.models:

  an object of class `candidateModels` containing a list of candidate
  models for stacking. See
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  for details.

- n.samples:

  number of samples to be drawn from the posterior distribution.

- loopd.controls:

  a list with details on how leave-one-out predictive densities (LOO-PD)
  are to be calculated. Valid tags include `method`, `CV.K` and `nMC`.
  The tag `method` can be either `'exact'` or `'CV'`. If sample size is
  more than 100, then the default is `'CV'` with `CV.K` equal to its
  default value 10 (Gelman *et al.* 2024). The tag `nMC` decides how
  many Monte Carlo samples will be used to evaluate the leave-one-out
  predictive densities, which must be at least 500 (default).

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

An object of class `stvcGLMstack`, which is a list including the
following tags -

- `samples`:

  a list of length equal to total number of candidate models with each
  entry corresponding to a list of length 3, containing posterior
  samples of fixed effects (`beta`), spatial effects (`z`), and fine
  scale variation `xi` for that model.

- `loopd`:

  a list of length equal to total number of candidate models with each
  entry containing leave-one-out predictive densities under that
  particular model.

- `n.models`:

  number of candidate models that are fit.

- `candidate.models`:

  a list of length `n_model` rows with each entry containing details of
  the model parameters.

- `stacking.weights`:

  a numeric vector of length equal to the number of candidate models
  storing the optimal stacking weights.

- `run.time`:

  a `proc_time` object with runtime details.

- `solver.status`:

  solver status as returned by the optimization routine.

This object can be further used to recover posterior samples of the
scale parameters in the model, and subsequrently, to make predictions at
new locations or times using the function
[`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md).

## Examples

``` r
# \donttest{
set.seed(1234)
data("sim_stvcPoisson")
dat <- sim_stvcPoisson[1:100, ]

# create list of candidate models (multivariate)
mod.list2 <- candidateModels(list(phi_s = list(2, 3),
                                  phi_t = list(1, 2),
                                  boundary = c(0.5, 0.75)), "cartesian")

# fit a spatial-temporal varying coefficient model using predictive stacking
mod1 <- stvcGLMstack(y ~ x1 + (x1), data = dat, family = "poisson",
                     sp_coords = as.matrix(dat[, c("s1", "s2")]),
                     time_coords = as.matrix(dat[, "t_coords"]),
                     cor.fn = "gneiting-decay",
                     process.type = "multivariate",
                     candidate.models = mod.list2,
                     loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
                     n.samples = 500)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.8.1 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.051s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -246.726
#> ℹ Compile time: 0.051s
#> ℹ Solver time: 0.008s
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi_s | phi_t | boundary | weight |
#> +---------+-------+-------+----------+--------+
#> | Model 1 |      2|      1|      0.50| 0.000  |
#> | Model 2 |      3|      1|      0.50| 0.244  |
#> | Model 3 |      2|      2|      0.50| 0.000  |
#> | Model 4 |      3|      2|      0.50| 0.100  |
#> | Model 5 |      2|      1|      0.75| 0.000  |
#> | Model 6 |      3|      1|      0.75| 0.000  |
#> | Model 7 |      2|      2|      0.75| 0.000  |
#> | Model 8 |      3|      2|      0.75| 0.657  |
#> +---------+-------+-------+----------+--------+
#> 
# }
```
