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
data(simSpatial)
dat <- simSpatial[1:100, ]
cand.mod <- candidateModels(list(phi = c(3, 6), nu = c(0.5, 1),
                                 boundary = c(0.5)), "cartesian")

mod1 <- spGLMstack(y_pois ~ x1 + x2, data = dat, family = "poisson",
                   coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
                   candidate.models = cand.mod,
                   n.samples = 100,
                   loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
                   verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>           | phi | nu  | boundary | weight |
#> +---------+-----+-----+----------+--------+
#> | Model 1 |    3|  0.5|       0.5| 0      |
#> | Model 2 |    6|  0.5|       0.5| 0      |
#> | Model 3 |    3|  1.0|       0.5| 0      |
#> | Model 4 |    6|  1.0|       0.5| 1      |
#> +---------+-----+-----+----------+--------+
#> 

# Recover posterior samples of scale parameters
mod1.1 <- recoverGLMscale(mod1)

# sample from the stacked posterior distribution
post_samps <- stackedSampler(mod1.1)
```
