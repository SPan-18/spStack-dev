# Bayesian spatially varying coefficients generalized linear model

Fits a Bayesian generalized linear model with spatially varying
coefficients under fixed values of the spatial process parameters and
some auxiliary model parameters. The output contains posterior samples
of the fixed effects, the spatially varying coefficients and, if
required, finds leave-one-out predictive densities.

## Usage

``` r
svcGLMexact(
  formula,
  data = parent.frame(),
  family,
  coords,
  cor.fn,
  process.type,
  spParams,
  priors,
  boundary = 0.5,
  n.samples,
  loopd = FALSE,
  loopd.method = "exact",
  CV.K = 10,
  loopd.nMC = 500,
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
  typically the environment from which `svcGLMexact` is called.

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
  details.

- spParams:

  fixed values of the spatial process parameters, a list with tags `phi`
  (spatial decay) and, if `cor.fn = 'matern'`, `nu` (smoothness). If
  `process.type = 'independent'`, each is a vector of length \\r\\,
  otherwise a scalar.

- priors:

  (optional) a list with each tag corresponding to a hyperparameter name
  and containing hyperprior details. Valid tags include `V.beta`,
  `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
  and `nu.z` must be at least 2.1 (`nu.z`, the inverse-Wishart degrees
  of freedom if `process.type = 'multivariate'`, at least 1). If not
  supplied, uses defaults.

- boundary:

  Specifies the boundary adjustment parameter. Must be a real number
  between 0 and 1. Default is 0.5.

- n.samples:

  number of posterior samples to be generated.

- loopd:

  logical. If `loopd=TRUE`, returns leave-one-out predictive densities,
  using method as given by `loopd.method`. Default is `FALSE`.

- loopd.method:

  character. Ignored if `loopd=FALSE`. If `loopd=TRUE`, valid inputs are
  `'exact'`, `'CV'`. The option `'exact'` corresponds to exact
  leave-one-out predictive densities which requires computation almost
  equivalent to fitting the model \\n\\ times. The options `'CV'` is
  faster as it implements \\K\\-fold cross validation to find
  approximate leave-one-out predictive densities (Vehtari *et al.*
  2017).

- CV.K:

  An integer between 10 and 20. Considered only if `loopd.method='CV'`.
  Default is 10 (as recommended in Vehtari *et. al* 2017).

- loopd.nMC:

  Number of Monte Carlo samples to be used to evaluate leave-one-out
  predictive densities when `loopd.method` is set to either 'exact' or
  'CV'.

- verbose:

  logical. If `verbose = TRUE`, prints model description.

- ...:

  currently no additional argument.

## Value

An object of class `svcGLMexact`, which is a list with the following
tags -

- priors:

  details of the priors used, containing the values of the boundary
  adjustment parameter (`boundary`), the variance parameter of the
  fine-scale variation term (`sigmaSq.xi`) and others.

- samples:

  a list of length 3, containing posterior samples of fixed effects
  (`beta`), the spatially varying coefficients (`z`, an \\nr \times\\
  `n.samples` matrix whose rows \\(j-1)n + 1, \ldots, jn\\ correspond to
  the \\j\\-th varying coefficient) and the fine-scale variation term
  (`xi`).

- loopd:

  If `loopd=TRUE`, contains leave-one-out predictive densities.

- model.params:

  Values of the fixed parameters: `phi` and `nu` (`NA` for the
  exponential correlation function).

- diagnostics:

  a list of diagnostics. Element `numerical` is a data frame with one
  row (one row per process, if `process.type = 'independent'`) and
  columns `min.pivot` (the smallest relative Cholesky pivot of the \\n
  \times n\\ spatial correlation matrix; values below 1e-8 indicate a
  nearly singular correlation matrix), `min.cor` and `max.cor` (the
  correlations of the two farthest-apart and of the two closest
  locations; values of `min.cor` above 0.95 suggest an effective range
  far exceeding the extent of the data, values of `max.cor` below 0.05
  nearly uncorrelated locations). If `verbose = TRUE`, a "Diagnostics"
  section is printed when any threshold is crossed.

The return object might include additional data that can be used for
subsequent prediction and/or model fit evaluation.

## Details

With this function, we fit a Bayesian hierarchical spatially varying
coefficients generalized linear model by sampling exactly from the joint
posterior distribution utilizing the generalized conjugate multivariate
distribution theory (Bradley and Clinch 2024). Suppose \\\chi = (s_1,
\ldots, s_n)\\ denotes the \\n\\ spatial locations in \\\mathcal{S}
\subset \mathbb{R}^2\\ at which the response \\y\\ is observed, with
\\y(s)\\ endowed with a probability law from the natural exponential
family \$\$ y(s) \sim \mathrm{EF}(x(s)^\top \beta + \tilde{x}(s)^\top
z(s); b(s), \psi) \$\$ for some positive parameter \\b(s) \> 0\\ and
unit log partition function \\\psi\\. Here, \\\tilde{x}(s)\\ denotes the
\\r\\ covariates with spatially varying coefficients. The supported
families (`'poisson'`, `'binomial'` and `'binary'`) and the hierarchical
model are those of
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
with \$\$ z_j \mid \sigma^2\_{z_j} \sim N(0, \sigma^2\_{z_j} R(\chi;
\phi_j, \nu_j)), \quad \sigma^2\_{z_j} \sim \mathrm{IG}(\nu_z/2,
\nu_z/2), \quad j = 1, \ldots, r, \$\$ where \\R(\chi; \phi, \nu)\\ is
the exponential or Matern spatial correlation matrix with decay \\\phi\\
and smoothness \\\nu\\. If `process.type = 'independent'`, each process
has its own \\(\phi_j, \nu_j)\\; if
`process.type = 'independent.shared'`, they share one. If
`process.type = 'multivariate'`, \\z = (z_1, \ldots, z_r)\\ is a
multivariate Gaussian process with correlation matrix \\R(\chi; \phi,
\nu)\\ and an \\r \times r\\ cross-covariance matrix with an
inverse-Wishart prior. We fix the process parameters and the
hyperparameters \\V\_\beta\\, \\\nu\_\beta\\, \\\nu_z\\ and
\\\sigma^2\_\xi\\; their defaults are as in
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md).

## References

Bradley JR, Clinch M (2024). "Generating Independent Replicates Directly
from the Posterior Distribution for a Class of Spatial Hierarchical
Models." *Journal of Computational and Graphical Statistics*, **0**(0),
1-17.
[doi:10.1080/10618600.2024.2365728](https://doi.org/10.1080/10618600.2024.2365728)
.

Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian Inference for
Spatial-temporal Non-Gaussian Data Using Predictive Stacking." *Bayesian
Analysis*, **In Press**.
[doi:10.1214/25-BA1582](https://doi.org/10.1214/25-BA1582) .

Vehtari A, Gelman A, Gabry J (2017). "Practical Bayesian Model
Evaluation Using Leave-One-out Cross-Validation and WAIC." *Statistics
and Computing*, **27**(5), 1413-1432. ISSN 0960-3174.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)
.

## See also

[`svcGLMstack()`](https://span-18.github.io/spStack-dev/reference/svcGLMstack.md),
[`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md),
[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
data(simSpatial)
dat <- simSpatial[1:150, ]

mod1 <- svcGLMexact(y_svc_pois ~ x1 + x2 + (x1), data = dat,
                    family = "poisson",
                    coords = as.matrix(dat[, c("s1", "s2")]),
                    cor.fn = "matern", process.type = "independent",
                    spParams = list(phi = c(8, 4), nu = c(1.5, 1.5)),
                    n.samples = 100, verbose = FALSE)

# rows 1, ..., n of z hold the varying intercept, rows n+1, ..., 2n the
# varying slope of x1
n <- nrow(dat)
z_hat <- apply(mod1$samples$z, 1, median)
cor(z_hat[n + 1:n], dat$z_true)
#> [1] 0.8783756
```
