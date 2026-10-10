# Univariate Bayesian spatial generalized linear model

Fits a Bayesian spatial generalized linear model with fixed values of
spatial process parameters and some auxiliary model parameters. The
output contains posterior samples of the fixed effects, spatial random
effects and, if required, finds leave-one-out predictive densities.

## Usage

``` r
spGLMexact(
  formula,
  data = parent.frame(),
  family,
  coords,
  cor.fn,
  priors,
  spParams,
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

  a symbolic description of the regression model to be fit. See example
  below.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `spGLMexact` is called.

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

  (optional) a list with each tag corresponding to a hyperparameter name
  and containing hyperprior details. Valid tags include `V.beta`,
  `nu.beta`, `nu.z` and `sigmaSq.xi`. Values of `nu.beta` and `nu.z`
  must be at least 2.1. If not supplied, uses defaults.

- spParams:

  fixed values of spatial process parameters.

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
  `'exact'`, `'CV'` and `'PSIS'`. The option `'exact'` corresponds to
  exact leave-one-out predictive densities which requires computation
  almost equivalent to fitting the model \\n\\ times. The options `'CV'`
  and `'PSIS'` are faster and they implement \\K\\-fold cross validation
  and Pareto-smoothed importance sampling to find approximate
  leave-one-out predictive densities (Vehtari *et al.* 2017).

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

An object of class `spGLMexact`, which is a list with the following tags
-

- priors:

  details of the priors used, containing the values of the boundary
  adjustment parameter (`boundary`), the variance parameter of the
  fine-scale variation term (`simasq.xi`) and others.

- samples:

  a list of length 3, containing posterior samples of fixed effects
  (`beta`), spatial effects (`z`) and the fine-scale variation term
  (`xi`).

- loopd:

  If `loopd=TRUE`, contains leave-one-out predictive densities.

- model.params:

  Values of the fixed parameters that includes `phi` (spatial decay),
  `nu` (spatial smoothness; `NA` for the exponential correlation
  function).

- diagnostics:

  a list of diagnostics. Element `numerical` is a data frame with one
  row and columns `min.pivot` (the smallest relative Cholesky pivot of
  the \\n \times n\\ correlation matrix; values below 1e-8 indicate a
  nearly singular correlation matrix), `min.cor` and `max.cor` (the
  correlations of the two farthest-apart and of the two closest
  locations; values of `min.cor` above 0.95 suggest an effective range
  far exceeding the extent of the data, values of `max.cor` below 0.05
  nearly uncorrelated locations), obtained from quantities the fit
  computes anyway. If `verbose = TRUE`, a "Diagnostics" section is
  printed when any threshold is crossed.

The return object might include additional data that can be used for
subsequent prediction and/or model fit evaluation.

## Details

With this function, we fit a Bayesian hierarchical spatial generalized
linear model by sampling exactly from the joint posterior distribution
utilizing the generalized conjugate multivariate distribution theory
(Bradley and Clinch 2024). Suppose \\\chi = (s_1, \ldots, s_n)\\ denotes
the \\n\\ spatial locations the response \\y\\ is observed. Let \\y(s)\\
be the outcome at location \\s\\ endowed with a probability law from the
natural exponential family, which we denote by \$\$ y(s) \sim
\mathrm{EF}(x(s)^\top \beta + z(s); b, \psi) \$\$ for some positive
parameter \\b \> 0\\ and unit log partition function \\\psi\\. We
consider the following response models based on the input supplied to
the argument `family`.

- `'poisson'`:

  It considers point-referenced Poisson responses \\y(s) \sim
  \mathrm{Poisson}(e^{x(s)^\top \beta + z(s)})\\. Here, \\b = 1\\ and
  \\\psi(t) = e^t\\.

- `'binomial'`:

  It considers point-referenced binomial counts \\y(s) \sim
  \mathrm{Binomial}(m(s), \pi(s))\\ where, \\m(s)\\ denotes the total
  number of trials and probability of success \\\pi(s) =
  \mathrm{ilogit}(x(s)^\top \beta + z(s))\\ at location \\s\\. Here, \\b
  = m(s)\\ and \\\psi(t) = \log(1+e^t)\\.

- `'binary'`:

  It considers point-referenced binary data (0 or, 1) i.e., \\y(s) \sim
  \mathrm{Bernoulli}(\pi(s))\\, where probability of success \\\pi(s) =
  \mathrm{ilogit}(x(s)^\top \beta + z(s))\\ at location \\s\\. Here, \\b
  = 1\\ and \\\psi(t) = \log(1 + e^t)\\.

The hierarchical model is given as \$\$ \begin{aligned} y(s_i) &\mid
\beta, z, \xi \sim EF(x(s_i)^\top \beta + z(s_i) + \xi_i - \mu_i; b_i,
\psi_y), i = 1, \ldots, n\\ \xi &\mid \beta, z, \sigma^2\_\xi,
\alpha\_\epsilon \sim \mathrm{GCM_c}(\cdots),\\ \beta &\mid
\sigma^2\_\beta \sim N(0, \sigma^2\_\beta V\_\beta), \quad
\sigma^2\_\beta \sim \mathrm{IG}(\nu\_\beta/2, \nu\_\beta/2)\\ z &\mid
\sigma^2_z \sim N(0, \sigma^2_z R(\chi; \phi, \nu)), \quad \sigma^2_z
\sim \mathrm{IG}(\nu_z/2, \nu_z/2), \end{aligned} \$\$ where \\\mu =
(\mu_1, \ldots, \mu_n)^\top\\ denotes the discrepancy parameter. We fix
the spatial process parameters \\\phi\\ and \\\nu\\ and the
hyperparameters \\V\_\beta\\, \\\nu\_\beta\\, \\\nu_z\\ and
\\\sigma^2\_\xi\\. The term \\\xi\\ is known as the fine-scale variation
term which is given a conditional generalized conjugate multivariate
distribution as prior. For more details, see Pan *et al.* 2025. Default
values for \\V\_\beta\\, \\\nu\_\beta\\, \\\nu_z\\, \\\sigma^2\_\xi\\
are diagonal with each diagonal element 100, 2.1, 2.1 and 0.1
respectively.

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

[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)

## Author

Soumyakanti Pan <span18@ucla.edu>

## Examples

``` r
data(simSpatial)
dat <- simSpatial[1:100, ]
coords <- as.matrix(dat[, c("s1", "s2")])

# Example 1: spatial Poisson count data
mod1 <- spGLMexact(y_pois ~ x1 + x2, data = dat, family = "poisson",
                   coords = coords, cor.fn = "matern",
                   spParams = list(phi = 6, nu = 0.5),
                   n.samples = 100, verbose = TRUE)
#> ----------------------------------------
#>  Model description
#> ----------------------------------------
#> Model fit with 100 observations.
#> 
#> Family = poisson.
#> 
#> Number of covariates 3 (including intercept).
#> 
#> Using the matern spatial correlation function.
#> 
#> Priors:
#>  beta: Gaussian
#>  mu: 0.00    0.00    0.00    
#>  cov:
#>   100.00  0.00    0.00   
#>   0.00    100.00  0.00   
#>   0.00    0.00    100.00 
#> 
#>  sigmaSq.beta ~ IG(nu.beta/2, nu.beta/2)
#>  sigmaSq.z ~ IG(nu.z/2, nu.z/2)
#>  nu.beta = 2.10, nu.z = 2.10.
#>  sigmaSq.xi = 0.10.
#>  Boundary adjustment parameter = 0.50.
#> 
#> Spatial process parameters:
#>  phi = 6.00, and, nu = 0.50.
#> 
#> Number of posterior samples = 100.
#> ----------------------------------------
post_beta <- mod1$samples$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>             2.5%        50%      97.5%
#> [1,]  0.85321873  1.7792510  3.5396481
#> [2,] -0.84012204 -0.5900074 -0.3406932
#> [3,]  0.07424866  0.4209770  0.7270565

# Example 2: spatial binomial count data
mod2 <- spGLMexact(cbind(y_binom, n_trials) ~ x1 + x2, data = dat,
                   family = "binomial", coords = coords, cor.fn = "matern",
                   spParams = list(phi = 6, nu = 0.5),
                   n.samples = 100, verbose = FALSE)

# Example 3: spatial binary data
mod3 <- spGLMexact(y_bin ~ x1 + x2, data = dat, family = "binary",
                   coords = coords, cor.fn = "matern",
                   spParams = list(phi = 6, nu = 0.5),
                   n.samples = 100, verbose = FALSE)
```
