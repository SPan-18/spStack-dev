# Bayesian spatially varying coefficients linear model

Fits a Bayesian linear model with spatially varying coefficients for a
Gaussian response, with the spatial process parameters and the
noise-to-spatial variance ratios fixed to values supplied by the user.
The output contains exact posterior samples of the fixed effects, the
noise and process variances, the spatially varying coefficients and, if
required, leave-one-out predictive densities.

## Usage

``` r
svcLMexact(
  formula,
  data = parent.frame(),
  coords,
  cor.fn,
  process.type,
  spParams,
  noise_sp_ratio,
  priors = "flat",
  n.samples,
  loopd = FALSE,
  loopd.method = "exact",
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
  typically the environment from which `svcLMexact` is called.

- coords:

  an \\n \times 2\\ matrix of the observation coordinates in
  \\\mathbb{R}^2\\ (e.g., easting and northing).

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial dependence structure among the observations. Supported
  covariance model key words are: `'exponential'` and `'matern'`. See
  below for details.

- process.type:

  a quoted keyword specifying the model for the spatial processes of the
  varying coefficients. Supported keywords are `'independent'`,
  independent processes with their own process parameters and
  noise-to-spatial variance ratios, and `'independent.shared'`,
  independent processes that share common process parameters and a
  common noise-to-spatial variance ratio.

- spParams:

  fixed values of the spatial process parameters, a list with tags `phi`
  (spatial decay) and, if `cor.fn = 'matern'`, `nu` (smoothness). If
  `process.type = 'independent'`, each is a vector of length \\r\\,
  otherwise a scalar.

- noise_sp_ratio:

  noise-to-spatial variance ratio(s) \\\delta^2_j\\: a vector of length
  \\r\\ if `process.type = 'independent'`, otherwise a scalar. Default
  is 1.

- priors:

  either `"flat"` (default), which assigns the prior \\p(\beta,
  \sigma^2) \propto 1/\sigma^2\\, or a list with tags `beta.norm` (a
  list containing \\\mu\_\beta\\ and \\V\_\beta\\) and/or `sigma.sq.ig`
  (a vector containing \\a\_\sigma\\ and \\b\_\sigma\\). A component not
  supplied in the list receives its flat prior, \\p(\beta) \propto 1\\
  or \\p(\sigma^2) \propto 1/\sigma^2\\.

- n.samples:

  number of posterior samples to be generated.

- loopd:

  logical. If `loopd=TRUE`, returns leave-one-out predictive densities,
  using method as given by `loopd.method`. Default is `FALSE`.

- loopd.method:

  character. Ignored if `loopd=FALSE`. If `loopd=TRUE`, valid inputs are
  `'exact'` and `'PSIS'`. The option `'exact'` finds the exact
  leave-one-out predictive densities in closed form, at the cost of
  about one additional \\n \times n\\ triangular inversion. The option
  `'PSIS'` finds approximate leave-one-out predictive densities using
  Pareto-smoothed importance sampling (Vehtari *et al.* 2024); with many
  latent effects (\\nr\\), its Pareto \\k\\ diagnostics are often high
  and `'exact'` is recommended.

- verbose:

  logical. If `verbose = TRUE`, prints model description.

- ...:

  currently no additional argument.

## Value

An object of class `svcLMexact`, which is a list with the following tags
-

- samples:

  a list of length 4, containing posterior samples of fixed effects
  (`beta`, a \\p \times\\ `n.samples` matrix), the noise variance
  (`sigmaSq`), the process variances (`sigmaSq.z`, an \\r \times\\
  `n.samples` matrix if `process.type = 'independent'`, otherwise a
  vector), and the spatially varying coefficients (`z`, an \\nr \times\\
  `n.samples` matrix whose rows \\(j-1)n + 1, \ldots, jn\\ correspond to
  the \\j\\-th varying coefficient).

- loopd:

  If `loopd=TRUE`, contains leave-one-out predictive densities.

- model.params:

  Values of the fixed parameters: `phi`, `nu` (`NA` for the exponential
  correlation function) and `noise_sp_ratio`.

- diagnostics:

  a list of fit diagnostics, obtained from quantities the fit computes
  anyway. Element `numerical` is a data frame with one row (one row per
  process, if `process.type = 'independent'`) and columns `min.pivot`
  (the smallest relative Cholesky pivot of the spatial correlation
  matrix and of \\V_y\\; values below 1e-8 indicate a nearly singular
  matrix), `min.cor` and `max.cor` (the correlations of the two
  farthest-apart and of the two closest locations; values of `min.cor`
  above 0.95 suggest an effective range far exceeding the extent of the
  data, values of `max.cor` below 0.05 nearly uncorrelated locations).
  If `loopd.method = 'PSIS'`, element `pareto` is a list with the Pareto
  \\k\\ diagnostic values (`k`), the threshold above which they are
  unreliable (`threshold`) and the number of values above it (`n.high`).
  If `verbose = TRUE`, a "Diagnostics" section is printed when any
  threshold is crossed.

The return object might include additional data used for subsequent
prediction and/or model fit evaluation.

## Details

Suppose \\\chi = (s_1, \ldots, s_n)\\ denotes the \\n\\ spatial
locations in \\\mathcal{S} \subset \mathbb{R}^2\\ at which the response
\\y\\ is observed. With this function, we fit the conjugate Bayesian
hierarchical model \$\$ \begin{aligned} y(s) &= x(s)^\top \beta +
\tilde{x}(s)^\top z(s) + \epsilon(s), \quad \epsilon(s) \sim N(0,
\sigma^2),\\ z_j &\mid \sigma^2 \sim N(0, \sigma^2\_{z_j} R(\chi;
\phi_j, \nu_j)), \quad \sigma^2\_{z_j} = \sigma^2 / \delta^2_j, \quad j
= 1, \ldots, r,\\ \beta &\mid \sigma^2 \sim N(\mu\_\beta, \sigma^2
V\_\beta), \quad \sigma^2 \sim \mathrm{IG}(a\_\sigma, b\_\sigma),
\end{aligned} \$\$ where \\\tilde{x}(s)\\ denotes the \\r\\ covariates
with spatially varying coefficients, the processes \\z_1, \ldots, z_r\\
are independent, and \\R(\chi; \phi, \nu)\\ is the exponential or Matern
spatial correlation matrix with decay \\\phi\\ and smoothness \\\nu\\.
We fix the noise-to-spatial variance ratios \\\delta^2_j = \sigma^2 /
\sigma^2\_{z_j}\\, the process parameters \\\phi_j\\ and \\\nu_j\\, and
the hyperparameters \\\mu\_\beta\\, \\V\_\beta\\, \\a\_\sigma\\ and
\\b\_\sigma\\. If `process.type = 'independent'`, each process has its
own \\(\phi_j, \nu_j, \delta^2_j)\\; if
`process.type = 'independent.shared'`, they share one. If
`priors = "flat"`, we instead assign the prior \\p(\beta, \sigma^2)
\propto 1/\sigma^2\\.

The joint posterior distribution is available in closed form and is
sampled exactly by composition, \$\$ p(\sigma^2, \beta, z \mid y) =
p(\sigma^2 \mid y) \times p(\beta \mid \sigma^2, y) \times p(z \mid
\beta, \sigma^2, y), \$\$ where \\\sigma^2 \mid y\\ is inverse-gamma and
the other two are Gaussian. All of them depend on the data through the
\\n \times n\\ matrix \\V_y = I_n + \sum_j \delta_j^{-2} D_j R_j D_j\\,
with \\D_j = \mathrm{diag}(\tilde{x}\_j)\\. The \\nr\\-dimensional
vector \\z\\ is drawn by a prior draw followed by a kriging correction
(Matheron's rule; Bhattacharya, Chakraborty and Mallick 2016), which
needs only \\n \times n\\ Cholesky factorizations. The exact
leave-one-out predictive densities are obtained in closed form from the
same factorizations. The computations are those of
[`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md),
with a spatial correlation function in place of the spatial-temporal
one.

## References

Bhattacharya A, Chakraborty A, Mallick BK (2016). "Fast sampling with
Gaussian scale mixture priors in high-dimensional regression."
*Biometrika*, **103**(4), 985-991.
[doi:10.1093/biomet/asw042](https://doi.org/10.1093/biomet/asw042) .

Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto Smoothed
Importance Sampling." *Journal of Machine Learning Research*,
**25**(72), 1-58. URL <https://jmlr.org/papers/v25/19-556.html>.

## See also

[`svcLMstack()`](https://span-18.github.io/spStack-dev/reference/svcLMstack.md),
[`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md),
[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
[`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
data(simSpatial)
dat <- simSpatial[1:200, ]

# varying intercept and varying slope of x1, each with its own parameters
mod1 <- svcLMexact(y_svc ~ x1 + x2 + (x1), data = dat,
                   coords = as.matrix(dat[, c("s1", "s2")]),
                   cor.fn = "matern", process.type = "independent",
                   spParams = list(phi = c(8, 4), nu = c(1.5, 1.5)),
                   noise_sp_ratio = c(0.25, 0.25),
                   n.samples = 500, loopd = TRUE, verbose = FALSE)

# rows 1, ..., n of z hold the varying intercept, rows n+1, ..., 2n the
# varying slope of x1
n <- nrow(dat)
z_hat <- apply(mod1$samples$z, 1, median)
cor(z_hat[1:n], dat$z2_true)
#> [1] 0.9514149
cor(z_hat[n + 1:n], dat$z_true)
#> [1] 0.9692254
```
