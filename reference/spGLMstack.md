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
  params.list,
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

- params.list:

  a list containing candidate values of spatial process parameters for
  the `cor.fn` used, and, the boundary parameter.

- n.samples:

  number of posterior samples to be generated.

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
  leave-one-out predictive densities.

- `n.models`:

  number of candidate models that are fit.

- `candidate.models`:

  a matrix with `n_model` rows with each row containing details of the
  model parameters and its optimal weight.

- `stacking.weights`:

  a numeric vector of length equal to the number of candidate models
  storing the optimal stacking weights.

- `run.time`:

  a `proc_time` object with runtime details.

- `solver.status`:

  solver status as returned by the optimization routine.

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
data("simPoisson")
dat <- simPoisson[1:100,]
mod1 <- spGLMstack(y ~ x1, data = dat, family = "poisson",
                   coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
                  params.list = list(phi = c(3, 7, 10), nu = c(0.25, 0.5, 1.5),
                                     boundary = c(0.5, 0.6)),
                  n.samples = 1000,
                  loopd.controls = list(method = "CV", CV.K = 10, nMC = 1000),
                  parallel = TRUE, verbose = TRUE)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.8.1 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.038s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -157.72
#> ℹ Compile time: 0.038s
#> ℹ Solver time: 0.006s
#> 
#> STACKING WEIGHTS:
#> 
#>            | phi | nu   | boundary | weight |
#> +----------+-----+------+----------+--------+
#> | Model 1  |    3|  0.25|       0.5| 0.000  |
#> | Model 2  |    7|  0.25|       0.5| 0.000  |
#> | Model 3  |   10|  0.25|       0.5| 0.000  |
#> | Model 4  |    3|  0.50|       0.5| 0.000  |
#> | Model 5  |    7|  0.50|       0.5| 0.000  |
#> | Model 6  |   10|  0.50|       0.5| 0.000  |
#> | Model 7  |    3|  1.50|       0.5| 0.379  |
#> | Model 8  |    7|  1.50|       0.5| 0.000  |
#> | Model 9  |   10|  1.50|       0.5| 0.000  |
#> | Model 10 |    3|  0.25|       0.6| 0.000  |
#> | Model 11 |    7|  0.25|       0.6| 0.000  |
#> | Model 12 |   10|  0.25|       0.6| 0.000  |
#> | Model 13 |    3|  0.50|       0.6| 0.000  |
#> | Model 14 |    7|  0.50|       0.6| 0.000  |
#> | Model 15 |   10|  0.50|       0.6| 0.000  |
#> | Model 16 |    3|  1.50|       0.6| 0.000  |
#> | Model 17 |    7|  1.50|       0.6| 0.621  |
#> | Model 18 |   10|  1.50|       0.6| 0.000  |
#> +----------+-----+------+----------+--------+
#> 

# print(mod1$solver.status)
# print(mod1$run.time)

post_samps <- stackedSampler(mod1)
post_beta <- post_samps$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>                   2.5%        50%      97.5%
#> (Intercept)  0.1876210  2.1365567  4.8390871
#> x1          -0.6906675 -0.5686308 -0.4223164

post_z <- post_samps$z
post_z_summ <- t(apply(post_z, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))

z_combn <- data.frame(z = dat$z_true,
                      zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2],
                      zU = post_z_summ[, 3])

library(ggplot2)
plot_z <- ggplot(data = z_combn, aes(x = z)) +
 geom_errorbar(aes(ymin = zL, ymax = zU),
               width = 0.05, alpha = 0.15,
               color = "skyblue") +
 geom_point(aes(y = zM), size = 0.25,
            color = "darkblue", alpha = 0.5) +
 geom_abline(slope = 1, intercept = 0,
             color = "red", linetype = "solid") +
 xlab("True z") + ylab("Posterior of z") +
 theme_bw() +
 theme(panel.background = element_blank(),
       aspect.ratio = 1)
# }
```
