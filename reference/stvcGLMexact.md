# Bayesian spatially-temporally varying generalized linear model

Fits a Bayesian generalized linear model with spatially-temporally
varying coefficients under fixed values of spatial-temporal process
parameters and some auxiliary model parameters. The output contains
posterior samples of the fixed effects, spatial-temporal random effects
and, if required, finds leave-one-out predictive densities.

## Usage

``` r
stvcGLMexact(
  formula,
  data = parent.frame(),
  family,
  sp_coords,
  time_coords,
  cor.fn,
  process.type,
  sptParams,
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
  parenthesis are assigned spatially-temporally varying coefficients.
  See examples.

- data:

  an optional data frame containing the variables in the model. If not
  found in `data`, the variables are taken from `environment(formula)`,
  typically the environment from which `stvcGLMexact` is called.

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

- sptParams:

  fixed values of spatial-temporal process parameters in usually a list
  of length 2. If `cor.fn='gneiting-decay'`, then it is a list of length
  2 with tags `phi_s` and `phi_t`. If `process.type='independent'`, then
  `phi_s` and `phi_t` contain fixed values of the \\r\\ spatial-temporal
  processes, otherwise they will contain scalars. See examples below.

- priors:

  (optional) a list with each tag corresponding to a hyperparameter name
  and containing hyperprior details. Valid tags include `V.beta`,
  `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
  and `nu.z` must be at least 2.1. If not supplied, uses defaults.

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

An object of class `stvcGLMexact`, which is a list with the following
tags -

- priors:

  details of the priors used, containing the values of the boundary
  adjustment parameter (`boundary`), the variance parameter of the
  fine-scale variation term (`simasq.xi`) and others.

- samples:

  a list of length 3, containing posterior samples of fixed effects
  (`beta`), spatial-temporal effects (`z`) and the fine-scale variation
  term (`xi`). The element with tag `z` will again be a list of length
  \\r\\, each containing posterior samples of the spatial-temporal
  random effects corresponding to each varying coefficient.

- loopd:

  If `loopd=TRUE`, contains leave-one-out predictive densities.

- model.params:

  Values of the fixed parameters that includes `phi_s` (spatial decay),
  `phi_t` (temporal smoothness).

The return object might include additional data that can be used for
subsequent prediction and/or model fit evaluation.

## Details

With this function, we fit a Bayesian hierarchical spatially-temporally
varying generalized linear model by sampling exactly from the joint
posterior distribution utilizing the generalized conjugate multivariate
distribution theory (Bradley and Clinch 2024). Suppose \\\chi = (\ell_1,
\ldots, \ell_n)\\ denotes the \\n\\ spatial-temporal co-ordinates in
\\\mathcal{L} = \mathcal{S} \times \mathcal{T}\\, the response \\y\\ is
observed. Let \\y(\ell)\\ be the outcome at the co-ordinate \\\ell\\
endowed with a probability law from the natural exponential family,
which we denote by \$\$ y(\ell) \sim \mathrm{EF}(x(\ell)^\top \beta +
\tilde{x}(\ell)^\top z(\ell); b(\ell), \psi) \$\$ for some positive
parameter \\b(\ell) \> 0\\ and unit log partition function \\\psi\\.
Here, \\\tilde{x}(\ell)\\ denotes covariates with spatially-temporally
varying coefficients We consider the following response models based on
the input supplied to the argument `family`.

- `'poisson'`:

  It considers point-referenced Poisson responses \\y(\ell) \sim
  \mathrm{Poisson}(e^{x(\ell)^\top \beta + \tilde{x}(\ell)^\top
  z(\ell)})\\. Here, \\b(\ell) = 1\\ and \\\psi(t) = e^t\\.

- `'binomial'`:

  It considers point-referenced binomial counts \\y(\ell) \sim
  \mathrm{Binomial}(m(\ell), \pi(\ell))\\ where, \\m(\ell)\\ denotes the
  total number of trials and probability of success \\\pi(\ell) =
  \mathrm{ilogit}(x(\ell)^\top \beta + \tilde{x}(\ell)^\top z(\ell))\\
  at spatial-temporal co-ordinate \\\ell\\. Here, \\b = m(\ell)\\ and
  \\\psi(t) = \log(1+e^t)\\.

- `'binary'`:

  It considers point-referenced binary data (0 or, 1) i.e., \\y(\ell)
  \sim \mathrm{Bernoulli}(\pi(\ell))\\, where probability of success
  \\\pi(\ell) = \mathrm{ilogit}(x(\ell)^\top \beta +
  \tilde{x}(\ell)^\top z(\ell))\\ at spatial-temporal co-ordinate
  \\\ell\\. Here, \\b(\ell) = 1\\ and \\\psi(t) = \log(1 + e^t)\\.

The hierarchical model is given as \$\$ \begin{aligned} y(\ell_i) &\mid
\beta, z, \xi \sim EF(x(\ell_i)^\top \beta + \tilde{x}(\ell_i)^\top
z(s_i) + \xi_i - \mu_i; b_i, \psi_y), i = 1, \ldots, n\\ \xi &\mid
\beta, z, \sigma^2\_\xi, \alpha\_\epsilon \sim \mathrm{GCM_c}(\cdots),\\
\beta &\mid \sigma^2\_\beta \sim N(0, \sigma^2\_\beta V\_\beta), \quad
\sigma^2\_\beta \sim \mathrm{IG}(\nu\_\beta/2, \nu\_\beta/2)\\ z_j &\mid
\sigma^2\_{z_j} \sim N(0, \sigma^2\_{z_j} R(\chi; \phi_s, \phi_t)),
\quad \sigma^2\_{z_j} \sim \mathrm{IG}(\nu_z/2, \nu_z/2), j = 1, \ldots,
r \end{aligned} \$\$ where \\\mu = (\mu_1, \ldots, \mu_n)^\top\\ denotes
the discrepancy parameter. We fix the spatial-temporal process
parameters \\\phi_s\\ and \\\phi_t\\ and the hyperparameters
\\V\_\beta\\, \\\nu\_\beta\\, \\\nu_z\\ and \\\sigma^2\_\xi\\. The term
\\\xi\\ is known as the fine-scale variation term which is given a
conditional generalized conjugate multivariate distribution as prior.
For more details, see Pan *et al.* 2025. Default values for
\\V\_\beta\\, \\\nu\_\beta\\, \\\nu_z\\, \\\sigma^2\_\xi\\ are diagonal
with each diagonal element 100, 2.1, 2.1 and 0.1 respectively.

## References

Bradley JR, Clinch M (2024). "Generating Independent Replicates Directly
from the Posterior Distribution for a Class of Spatial Hierarchical
Models." *Journal of Computational and Graphical Statistics*, **0**(0),
1-17.
[doi:10.1080/10618600.2024.2365728](https://doi.org/10.1080/10618600.2024.2365728)
.

T. Gneiting and P. Guttorp (2010). "Continuous-parameter spatio-temporal
processes." In *A.E. Gelfand, P.J. Diggle, M. Fuentes, and P Guttorp,
editors, Handbook of Spatial Statistics*, Chapman & Hall CRC Handbooks
of Modern Statistical Methods, p 427–436. Taylor and Francis.

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

[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)

## Author

Soumyakanti Pan <span18@ucla.edu>

## Examples

``` r
data("sim_stvcPoisson")
dat <- sim_stvcPoisson[1:100, ]

# Fit a spatial-temporal varying coefficient Poisson GLM
mod1 <- stvcGLMexact(y ~ x1 + (x1), data = dat, family = "poisson",
                     sp_coords = as.matrix(dat[, c("s1", "s2")]),
                     time_coords = as.matrix(dat[, "t_coords"]),
                     cor.fn = "gneiting-decay",
                     process.type = "multivariate",
                     sptParams = list(phi_s = 1, phi_t = 1),
                     verbose = FALSE, n.samples = 100)
```
