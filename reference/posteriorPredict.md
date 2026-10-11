# Prediction of latent process at new spatial or temporal locations

A function to sample from the posterior predictive distribution of the
latent spatial or spatial-temporal process.

## Usage

``` r
posteriorPredict(mod_out, coords_new, covars_new, joint = FALSE, nBinom_new)
```

## Arguments

- mod_out:

  an object returned by any model fit under fixed hyperparameters or
  using predictive stacking, i.e.,
  [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md),
  [`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md)
  or
  [`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md).

- coords_new:

  new spatial coordinates (an \\n\_{new} \times 2\\ matrix) or, for the
  spatial-temporal models, a list with tags `sp` (an \\n\_{new} \times
  2\\ matrix) and `time` (a vector or one-column matrix) at which the
  latent process, the mean, and the response is to be predicted.

- covars_new:

  new covariates at the new coordinates: a matrix or, for the
  spatial-temporal models, a list with tags `fixed` (covariates with
  fixed effects) and `vc` (covariates with spatially-temporally varying
  coefficients). See examples for the structure of this list.

- joint:

  a logical value indicating whether to return the joint posterior
  predictive samples of the latent process at the new locations or
  times. Defaults to `FALSE`.

- nBinom_new:

  a vector of the number of trials for each new prediction location or
  time. Only required if the model family is `"binomial"`. Defaults to a
  vector of ones, indicating one trial for each new prediction.

## Value

A modified object with the class name preceeded by the identifier `pp`
separated by a `.`. For example, if input is of class `spLMstack`, then
the output of this prediction function would be `pp.spLMstack`. The
entry with the tag `samples` is updated and will include samples from
the posterior predictive distribution of the latent process, the mean,
and the response at the new locations or times. An entry with the tag
`prediction` is added and contains the new coordinates and covariates,
and whether the joint posterior predictive samples were requested.

## See also

[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md),
[`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md),
[`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1234)
# training and test data sizes
n_train <- 100
n_pred <- 10

# Example 1: Spatial linear model
# split the data into training and prediction sets
data(simSpatial)
dat_train <- simSpatial[1:n_train, ]
dat_pred <- simSpatial[n_train + 1:n_pred, ]

# fit a spatial linear model using predictive stacking
cand.mod <- candidateModels(list(phi = c(2, 3, 4), nu = c(1, 1.5, 1.75),
                                 noise_sp_ratio = c(0.1, 0.25)), "cartesian")
mod1 <- spLMstack(y_gauss ~ x1 + x2, data = dat_train,
                  coords = as.matrix(dat_train[, c("s1", "s2")]),
                  cor.fn = "matern",
                  candidate.models = cand.mod,
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>            | phi | nu   | noise_sp_ratio | weight |
#> +----------+-----+------+----------------+--------+
#> | Model 1  |    2|  1.00|            0.10| 0.000  |
#> | Model 2  |    3|  1.00|            0.10| 0.000  |
#> | Model 3  |    4|  1.00|            0.10| 0.000  |
#> | Model 4  |    2|  1.50|            0.10| 0.899  |
#> | Model 5  |    3|  1.50|            0.10| 0.000  |
#> | Model 6  |    4|  1.50|            0.10| 0.000  |
#> | Model 7  |    2|  1.75|            0.10| 0.101  |
#> | Model 8  |    3|  1.75|            0.10| 0.000  |
#> | Model 9  |    4|  1.75|            0.10| 0.000  |
#> | Model 10 |    2|  1.00|            0.25| 0.000  |
#> | Model 11 |    3|  1.00|            0.25| 0.000  |
#> | Model 12 |    4|  1.00|            0.25| 0.000  |
#> | Model 13 |    2|  1.50|            0.25| 0.000  |
#> | Model 14 |    3|  1.50|            0.25| 0.000  |
#> | Model 15 |    4|  1.50|            0.25| 0.000  |
#> | Model 16 |    2|  1.75|            0.25| 0.000  |
#> | Model 17 |    3|  1.75|            0.25| 0.000  |
#> | Model 18 |    4|  1.75|            0.25| 0.000  |
#> +----------+-----+------+----------------+--------+
#> 

# prepare new coordinates and covariates for prediction
sp_pred <- as.matrix(dat_pred[, c("s1", "s2")])
X_new <- cbind(1, dat_pred$x1, dat_pred$x2)

# carry out posterior prediction
mod.pred <- posteriorPredict(mod1, coords_new = sp_pred, covars_new = X_new,
                             joint = TRUE)

# sample from the stacked posterior and posterior predictive distribution
post_samps <- stackedSampler(mod.pred)

# compare the predicted spatial effects with the truth
z_pred <- apply(post_samps$z.pred, 1, median)
cor(z_pred, dat_pred$z_true)
#> [1] 0.9932671
```
