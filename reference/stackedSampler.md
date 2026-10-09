# Sample from the stacked posterior distribution

A helper function to sample from the stacked posterior distribution to
obtain final posterior samples that can be used for subsequent analysis.
This function applies on outputs of functions
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)
and
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md).

## Usage

``` r
stackedSampler(mod_out, n.samples)
```

## Arguments

- mod_out:

  an object that is an output of a model fit or a prediction task, i.e.,
  the class should be either `spLMstack`, 'pp.spLMstack', `spGLMstack`,
  `pp.spGLMstack`, `stvcGLMexact`, or `pp.stvcGLMexact`.

- n.samples:

  (optional) If missing, inherits the number of posterior samples from
  the original output. Otherwise, it specifies number of posterior
  samples to draw from the stacked posterior. If it exceeds the number
  of posterior draws used in the original function, then a message is
  thrown and the samples are obtained by resampling. We recommended
  running the original model fit/prediction with enough samples.

## Value

An object of class `stacked_posterior`, which is a list that includes
the following tags -

- beta:

  samples of the fixed effect from the stacked joint posterior.

- z:

  samples of the spatial random effects from the stacked joint
  posterior.

The list may also include other scale parameters corresponding to the
model.

## Details

After obtaining the optimal stacking weights \\\hat{w}\_1, \ldots,
\hat{w}\_G\\, posterior inference of quantities of interest subsequently
proceed from the *stacked* posterior, \$\$ \tilde{p}(\cdot \mid y) =
\sum\_{g = 1}^G \hat{w}\_g p(\cdot \mid y, M_g), \$\$ where
\\\mathcal{M} = \\M_1, \ldots, M_g\\\\ is the collection of candidate
models.

## See also

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
cand.mod <- candidateModels(list(phi = c(1.5, 3),
                                 nu = c(0.5, 1),
                                 noise_sp_ratio = c(1)), "cartesian")

mod1 <- spLMstack(y ~ x1, data = dat,
                  coords = as.matrix(dat[, c("s1", "s2")]),
                  cor.fn = "matern",
                  candidate.models = cand.mod,
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> --------------------------------------------------
#> Solver diagnostics:
#> Installed solvers: CLARABEL, SCS, OSQP, HIGHS
#> Requested solver: DEFAULT (CLARABEL -> ECOS -> SCS)
#> Solver search order: CLARABEL -> SCS
#> --------------------------------------------------
#> ────────────────────────────────── CVXR v1.9.2 ─────────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (DCP)
#> ℹ Compilation: "CLARABEL" via CVXR::FlipObjective -> CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::Clarabel_Solver
#> ℹ Compile time: 0.032s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: -54.6976
#> ℹ Compile time: 0.032s
#> ℹ Solver time: 0.005s
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | noise_sp_ratio | weight |
#> +---------+-----+-----+----------------+--------+
#> | Model 1 |  1.5|  0.5|               1| 0.000  |
#> | Model 2 |  3.0|  0.5|               1| 0.333  |
#> | Model 3 |  1.5|  1.0|               1| 0.000  |
#> | Model 4 |  3.0|  1.0|               1| 0.667  |
#> +---------+-----+-----+----------------+--------+
#> 
print(mod1$solver)
#> [1] "CVXR:CLARABEL"
print(mod1$solver.status)
#> [1] "optimal"
print(mod1$run.time)
#>    user  system elapsed 
#>   0.179   0.362   0.140 

post_samps <- stackedSampler(mod1)
post_beta <- post_samps$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#>                 2.5%      50%    97.5%
#> (Intercept) 1.697664 2.395785 3.090683
#> x1          4.865601 4.979866 5.093619
```
