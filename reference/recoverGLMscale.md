# Recover posterior samples of scale parameters of spatial/spatial-temporal generalized linear models

A function to recover posterior samples of scale parameters that were
marginalized out during model fit. This is only applicable for spatial
or, spatial-temporal generalized linear models. This function applies on
outputs of functions that fits a spatial/spatial-temporal generalized
linear model, such as
[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
and
[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md).

## Usage

``` r
recoverGLMscale(mod_out)
```

## Arguments

- mod_out:

  an object returned by a fitting a spatial or spatial-temporal GLM.

## Value

An object of the same class as input, and updates the list tagged
`samples` with the posterior samples of the scale parameters. The new
tags are `sigmasq.beta` and `z.scale`.

## See also

[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
data("simPoisson")
dat <- simPoisson[1:100, ]
mod1 <- spGLMstack(y ~ x1, data = dat, family = "poisson",
                   coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
                   params.list = list(phi = c(3, 5, 7), nu = c(0.5, 1.5),
                                      boundary = c(0.5)),
                   n.samples = 100,
                   loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
                   verbose = TRUE)
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
#> ✔ Optimal value: -154.519
#> ℹ Compile time: 0.051s
#> ℹ Solver time: 0.004s
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | boundary | weight |
#> +---------+-----+-----+----------+--------+
#> | Model 1 |    3|  0.5|       0.5| 0.000  |
#> | Model 2 |    5|  0.5|       0.5| 0.000  |
#> | Model 3 |    7|  0.5|       0.5| 0.000  |
#> | Model 4 |    3|  1.5|       0.5| 0.398  |
#> | Model 5 |    5|  1.5|       0.5| 0.082  |
#> | Model 6 |    7|  1.5|       0.5| 0.520  |
#> +---------+-----+-----+----------+--------+
#> 

# Recover posterior samples of scale parameters
mod1.1 <- recoverGLMscale(mod1)

# sample from the stacked posterior distribution
post_samps <- stackedSampler(mod1.1)
```
