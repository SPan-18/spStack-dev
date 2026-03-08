# Optimal stacking weights

Obtains optimal stacking weights given leave-one-out predictive
densities for each candidate model.

## Usage

``` r
get_stacking_weights(log_loopd, solver = NULL, verbose = TRUE)
```

## Arguments

- log_loopd:

  an \\n \times M\\ matrix with \\i\\-th row containing the
  leave-one-out predictive densities for the \\i\\-th data point for the
  \\M\\ candidate models.

- solver:

  specifies the solver to use for obtaining optimal weights. Default is
  `"CLARABEL"`. Internally calls
  [`CVXR::psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.html).

- verbose:

  if `TRUE`, prints output of optimization routine.

## Value

A list with elements:

- `weights`:

  optimal stacking weights as a numeric vector of length \\M\\

- `status`:

  solver status, returns `"optimal"` if solver succeeded.

- `solver`:

  name of the solver used.

## References

Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to Average
Bayesian Predictive Distributions (with Discussion)." *Bayesian
Analysis*, **13**(3), 917-1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091) .

## See also

[`CVXR::psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.html),
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
data(simGaussian)
dat <- simGaussian[1:100, ]

mod1 <- spLMstack(y ~ x1, data = dat,
                  coords = as.matrix(dat[, c("s1", "s2")]),
                  cor.fn = "matern",
                  params.list = list(phi = c(1.5, 3),
                                     nu = c(0.5, 1),
                                     noise_sp_ratio = c(1)),
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.8.1 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.824s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -50.8551
#> ℹ Compile time: 0.824s
#> ℹ Solver time: 0.025s
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | noise_sp_ratio | weight |
#> +---------+-----+-----+----------------+--------+
#> | Model 1 |  1.5|  0.5|               1| 0.000  |
#> | Model 2 |  3.0|  0.5|               1| 0.284  |
#> | Model 3 |  1.5|  1.0|               1| 0.000  |
#> | Model 4 |  3.0|  1.0|               1| 0.716  |
#> +---------+-----+-----+----------------+--------+
#> 

loopd_mat <- do.call('cbind', mod1$loopd)
w_hat <- get_stacking_weights(loopd_mat)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.8.1 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.214s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -50.8551
#> ℹ Compile time: 0.214s
#> ℹ Solver time: 0.057s
print(round(w_hat$weights, 4))
#> [1] 0.0000 0.2845 0.0000 0.7155
print(w_hat$solver)
#> [1] "CVXR:CLARABEL"
print(w_hat$status)
#> [1] "optimal"
```
