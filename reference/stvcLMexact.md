# Bayesian spatially-temporally varying coefficients linear model

Fits a Bayesian linear model with spatially-temporally varying
coefficients for a Gaussian response, with the spatial-temporal process
parameters and the noise-to-spatial variance ratios fixed to values
supplied by the user. The output contains exact posterior samples of the
fixed effects, the noise and process variances, the spatial-temporal
random effects and, if required, leave-one-out predictive densities.

## Usage

``` r
stvcLMexact(
  formula,
  data = parent.frame(),
  sp_coords,
  time_coords,
  cor.fn,
  process.type,
  sptParams,
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
  parenthesis are assigned spatially-temporally varying coefficients.
  See examples.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `stvcLMexact` is called.

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
  and Guttorp 2010).

- process.type:

  a quoted keyword specifying the model for the spatial-temporal
  processes of the varying coefficients. Supported keywords are
  `'independent'`, independent processes with their own process
  parameters and noise-to-spatial variance ratios, and
  `'independent.shared'`, independent processes that share common
  process parameters and a common noise-to-spatial variance ratio.

- sptParams:

  fixed values of the spatial-temporal process parameters, a list with
  tags `phi_s` and `phi_t`. If `process.type = 'independent'`, each is a
  vector of length \\r\\, otherwise a scalar.

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

An object of class `stvcLMexact`, which is a list with the following
tags -

- samples:

  a list of length 4, containing posterior samples of fixed effects
  (`beta`, a \\p \times\\ `n.samples` matrix), the noise variance
  (`sigmaSq`), the process variances (`sigmaSq.z`, an \\r \times\\
  `n.samples` matrix if `process.type = 'independent'`, otherwise a
  vector), and the spatial-temporal effects (`z`, an \\nr \times\\
  `n.samples` matrix whose rows \\(j-1)n + 1, \ldots, jn\\ correspond to
  the \\j\\-th varying coefficient).

- loopd:

  If `loopd=TRUE`, contains leave-one-out predictive densities.

- model.params:

  Values of the fixed parameters: `phi_s`, `phi_t` and `noise_sp_ratio`.

- diagnostics:

  a list of fit diagnostics, obtained from quantities the fit computes
  anyway. Element `numerical` is a data frame with one row (one row per
  process, if `process.type = 'independent'`) and columns `min.pivot`
  (the smallest relative Cholesky pivot of the spatial-temporal
  correlation matrix and of \\V_y\\; values below 1e-8 indicate a nearly
  singular matrix), `min.cor` and `max.cor` (the correlations of the two
  farthest-apart and of the two closest space-time locations; values of
  `min.cor` above 0.95 suggest an effective range far exceeding the
  extent of the data, values of `max.cor` below 0.05 nearly uncorrelated
  space-time locations). If `loopd.method = 'PSIS'`, element `pareto` is
  a list with the Pareto \\k\\ diagnostic values (`k`), the threshold
  above which they are unreliable (`threshold`) and the number of values
  above it (`n.high`). If `verbose = TRUE`, a "Diagnostics" section is
  printed when any threshold is crossed.

The return object might include additional data used for subsequent
prediction and/or model fit evaluation.

## Details

Suppose \\\chi = (\ell_1, \ldots, \ell_n)\\ denotes the \\n\\
spatial-temporal co-ordinates in \\\mathcal{L} = \mathcal{S} \times
\mathcal{T}\\ at which the response \\y\\ is observed. With this
function, we fit the conjugate Bayesian hierarchical model \$\$
\begin{aligned} y(\ell) &= x(\ell)^\top \beta + \tilde{x}(\ell)^\top
z(\ell) + \epsilon(\ell), \quad \epsilon(\ell) \sim N(0, \sigma^2),\\
z_j &\mid \sigma^2 \sim N(0, \sigma^2\_{z_j} R(\chi; \phi\_{s,j},
\phi\_{t,j})), \quad \sigma^2\_{z_j} = \sigma^2 / \delta^2_j, \quad j =
1, \ldots, r,\\ \beta &\mid \sigma^2 \sim N(\mu\_\beta, \sigma^2
V\_\beta), \quad \sigma^2 \sim \mathrm{IG}(a\_\sigma, b\_\sigma),
\end{aligned} \$\$ where \\\tilde{x}(\ell)\\ denotes the \\r\\
covariates with spatially-temporally varying coefficients, the processes
\\z_1, \ldots, z_r\\ are independent, and \\R(\chi; \phi_s, \phi_t)\\ is
the spatial-temporal correlation matrix of the Gneiting (2002) family.
We fix the noise-to-spatial variance ratios \\\delta^2_j = \sigma^2 /
\sigma^2\_{z_j}\\, the process parameters \\\phi\_{s,j}\\ and
\\\phi\_{t,j}\\, and the hyperparameters \\\mu\_\beta\\, \\V\_\beta\\,
\\a\_\sigma\\ and \\b\_\sigma\\. If `process.type = 'independent'`, each
process has its own \\(\phi\_{s,j}, \phi\_{t,j}, \delta^2_j)\\; if
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
needs only \\n \times n\\ Cholesky factorizations. Posterior samples of
the process variances are obtained as \\\sigma^2\_{z_j} = \sigma^2 /
\delta^2_j\\. The exact leave-one-out predictive densities are obtained
in closed form from the same factorizations.

## References

Bhattacharya A, Chakraborty A, Mallick BK (2016). "Fast sampling with
Gaussian scale mixture priors in high-dimensional regression."
*Biometrika*, **103**(4), 985-991.
[doi:10.1093/biomet/asw042](https://doi.org/10.1093/biomet/asw042) .

Gneiting T (2002). "Nonseparable, Stationary Covariance Functions for
Space-Time Data." *Journal of the American Statistical Association*,
**97**(458), 590-600.
[doi:10.1198/016214502760047113](https://doi.org/10.1198/016214502760047113)
.

T. Gneiting and P. Guttorp (2010). "Continuous-parameter spatio-temporal
processes." In *A.E. Gelfand, P.J. Diggle, M. Fuentes, and P Guttorp,
editors, Handbook of Spatial Statistics*, Chapman & Hall CRC Handbooks
of Modern Statistical Methods, p 427-436. Taylor and Francis.

Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto Smoothed
Importance Sampling." *Journal of Machine Learning Research*,
**25**(72), 1-58. URL <https://jmlr.org/papers/v25/19-556.html>.

## See also

[`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md),
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
data(simSpaceTime)
dat <- simSpaceTime[1:100, ]

mod1 <- stvcLMexact(y_gauss ~ x1 + x2 + (x1), data = dat,
                    sp_coords = as.matrix(dat[, c("s1", "s2")]),
                    time_coords = as.matrix(dat[, "t_coords"]),
                    cor.fn = "gneiting-decay",
                    process.type = "independent",
                    sptParams = list(phi_s = c(3, 6), phi_t = c(4, 2)),
                    noise_sp_ratio = c(0.5, 1),
                    n.samples = 500, loopd = TRUE, verbose = FALSE)

# rows 1, ..., n of z hold the varying intercept, rows n+1, ..., 2n the
# varying slope of x1
n <- nrow(dat)
z_hat <- apply(mod1$samples$z, 1, median)
cor(z_hat[1:n], dat$z1_true)
#> [1] 0.830181
```
