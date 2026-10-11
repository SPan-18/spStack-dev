# spStack: Bayesian Geostatistics Using Predictive Stacking

spStack implements Bayesian inference for a rich class of spatial and
spatial-temporal geostatistical models using stacking of predictive
densities. Besides delivering competitive predictive performance as
compared to traditional fully Bayesian inference using MCMC, predictive
stacking is embarrassingly parallel, and hence, fast. This package, to
the best of our knowledge, is the first to implement stacking for
Bayesian analysis of spatial and spatial-temporal data. Technical
details surrounding the methodology can be found in the articles Zhang
et al. ([2025](#ref-zhang2024stacking)) and Pan et al.
([2025](#ref-pan2024stacking)).

``` r

set.seed(1729)
```

## Bayesian Gaussian spatial regression models

Here is a quick example using the Gaussian response `y_gauss` of the
lazyloaded synthetic data `simSpatial`, observed at 500 locations in the
unit square with two covariates `x1` and `x2` and a spatial effect
`z_true` with a ring pattern (see
[`?simSpatial`](https://span-18.github.io/spStack-dev/reference/simSpatial.md)).

**Step 1.** Load the library spStack and prepare the data. Here we split
the data into `dat_train` and `dat_pred` - we train our model on
`dat_train` and subsequently carry out posterior predictive inference at
the locations in `dat_pred`.

``` r

library(spStack)

# training and test data sizes
n_train <- 200
n_pred <- 50

data("simSpatial")
dat_train <- simSpatial[1:n_train, ]
dat_pred <- simSpatial[n_train + 1:n_pred, ]
```

**Step 2.** Run the function
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md) -
define the model with a formula, input the spatial coordinates as a
matrix, specify the correlation function, and accordingly provide
candidate models through `candidate.models`. The argument `loopd.method`
can be used to specify the method used for calculation of leave-one-out
predictive densities, existing parallelization plan can be used if
`parallel` is set `TRUE`, and `solver` argument specifies the solver
used to carry out the optimization routine to get stacking weights.

``` r

cand.mod <- candidateModels(list(phi = c(2, 4, 6, 8),
                                 nu = c(0.5, 1.5, 2.5),
                                 noise_sp_ratio = c(0.1, 0.5)),
                            "cartesian")

mod1 <- spLMstack(y_gauss ~ x1 + x2, data = dat_train,
                  coords = as.matrix(dat_train[, c("s1", "s2")]),
                  cor.fn = "matern",
                  candidate.models = cand.mod,
                  n.samples = 1000, loopd.method = "exact",
                  parallel = FALSE, verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>            | phi | nu  | noise_sp_ratio | weight |
#> +----------+-----+-----+----------------+--------+
#> | Model 1  |    2|  0.5|             0.1| 0.000  |
#> | Model 2  |    4|  0.5|             0.1| 0.000  |
#> | Model 3  |    6|  0.5|             0.1| 0.000  |
#> | Model 4  |    8|  0.5|             0.1| 0.000  |
#> | Model 5  |    2|  1.5|             0.1| 0.000  |
#> | Model 6  |    4|  1.5|             0.1| 0.000  |
#> | Model 7  |    6|  1.5|             0.1| 0.000  |
#> | Model 8  |    8|  1.5|             0.1| 0.000  |
#> | Model 9  |    2|  2.5|             0.1| 0.411  |
#> | Model 10 |    4|  2.5|             0.1| 0.491  |
#> | Model 11 |    6|  2.5|             0.1| 0.000  |
#> | Model 12 |    8|  2.5|             0.1| 0.000  |
#> | Model 13 |    2|  0.5|             0.5| 0.000  |
#> | Model 14 |    4|  0.5|             0.5| 0.000  |
#> | Model 15 |    6|  0.5|             0.5| 0.000  |
#> | Model 16 |    8|  0.5|             0.5| 0.000  |
#> | Model 17 |    2|  1.5|             0.5| 0.000  |
#> | Model 18 |    4|  1.5|             0.5| 0.000  |
#> | Model 19 |    6|  1.5|             0.5| 0.000  |
#> | Model 20 |    8|  1.5|             0.5| 0.000  |
#> | Model 21 |    2|  2.5|             0.5| 0.000  |
#> | Model 22 |    4|  2.5|             0.5| 0.097  |
#> | Model 23 |    6|  2.5|             0.5| 0.000  |
#> | Model 24 |    8|  2.5|             0.5| 0.000  |
#> +----------+-----+-----+----------------+--------+
```

**Step 3.** Use the helper function
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md)
to sample from the stacked posterior distribution. These samples serve
as the final posterior samples corresponding to our target model.

``` r

post_samps <- stackedSampler(mod1)
```

The final output will be a tagged list with each entry containing
posterior samples of the corresponding parameter.

``` r

post_beta <- post_samps$beta
summary_beta <- t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod1$X.names
print(summary_beta)
#>                  2.5%        50%      97.5%
#> (Intercept) -4.841186 -0.6961954  1.5411549
#> x1           4.891173  4.9522738  5.0220084
#> x2          -1.057542 -0.9782068 -0.8994407
```

> **Note:** The following optional steps are only required if interested
> in posterior predictive inference.

**(Optional) Step 4.** Prepare the inputs for posterior predictive
inference. The new spatial coordinates at which we intend to carry out
prediction is given by `sp_pred` and the value of the covariates at
these new locations are given by `X_new`.

``` r

sp_pred <- as.matrix(dat_pred[, c("s1", "s2")])
X_new <- cbind(1, dat_pred$x1, dat_pred$x2)
```

**(Optional) Step 5.** Pass the output obtained by running
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)
through the function
[`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
along with the new coordinates and covariates.

``` r

mod.pred <- posteriorPredict(mod1,
                             coords_new = sp_pred,
                             covars_new = X_new,
                             joint = TRUE)
```

**(Optional) Step 6.** Once samples from the posterior predictive
distribution are obtained, once again run
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md)
to obtain samples from the *stacked* posterior predictive distribution.

``` r

postpred_samps <- stackedSampler(mod.pred)
```

Next, we analyze how well we predict the responses by plotting their
posterior predictive summaries against their corresponding true values.

``` r

postpred_y <- postpred_samps$y.pred
post_y_summ <- t(apply(postpred_y, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
y_combn <- data.frame(y = dat_pred$y_gauss, yL = post_y_summ[, 1],
                      yM = post_y_summ[, 2], yU = post_y_summ[, 3])
library(ggplot2)
ggplot(data = y_combn, aes(x = y)) +
  geom_errorbar(aes(ymin = yL, ymax = yU), color = "skyblue") +
  geom_point(aes(y = yM), color = "darkblue", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "solid") +
  xlab("True y") + ylab("Posterior predictive of y") + theme_bw() +
  theme(panel.background = element_blank(), aspect.ratio = 1)
```

![](spStack_files/figure-html/unnamed-chunk-10-1.png)

More details about this can found in the vignettes titled “Spatial
Models” and “Posterior Predictive Inference”.

## Bayesian non-Gaussian spatial regression models

The workflow for the spatial generalized linear models are similar. Here
is a quick example using the Poisson response `y_pois` of the same data
`simSpatial`.

**Step 1.** Prepare data by splitting into train and test sets.

``` r

# training and test data sizes
n_train <- 200
n_pred <- 50

dat_train <- simSpatial[1:n_train, ]
dat_pred <- simSpatial[n_train + 1:n_pred, ]
```

**Step 2.** Run the function
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
after specifying the `family`, and providing candidate models
constructed using
[`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md).
The `loopd.controls` option can be used to specify the method and
parameters used for calculation of leave-one-out predictive densities.
The input `list(method = "CV", CV.K = 10, nMC = 500)` corresponds to
$`K`$-fold cross-validation with $`K = 10`$ and using 500 Monte Carlo
samples for calculating each predictive density.

``` r

cand.mod <- candidateModels(list(phi = c(1, 2, 4, 8),
                                 nu = c(0.5, 1, 1.5, 2.5),
                                 boundary = c(0.5, 0.75)), "cartesian")

mod1 <- spGLMstack(y_pois ~ x1 + x2, data = dat_train, family = "poisson",
                   coords = as.matrix(dat_train[, c("s1", "s2")]), cor.fn = "matern",
                   candidate.models = cand.mod,
                   n.samples = 1000,
                   loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
                   verbose = TRUE)
#> 
#> STACKING WEIGHTS:
#> 
#>            | phi | nu  | boundary | weight |
#> +----------+-----+-----+----------+--------+
#> | Model 1  |    1|  0.5|      0.50| 0.000  |
#> | Model 2  |    2|  0.5|      0.50| 0.000  |
#> | Model 3  |    4|  0.5|      0.50| 0.000  |
#> | Model 4  |    8|  0.5|      0.50| 0.000  |
#> | Model 5  |    1|  1.0|      0.50| 0.000  |
#> | Model 6  |    2|  1.0|      0.50| 0.000  |
#> | Model 7  |    4|  1.0|      0.50| 0.000  |
#> | Model 8  |    8|  1.0|      0.50| 0.000  |
#> | Model 9  |    1|  1.5|      0.50| 0.000  |
#> | Model 10 |    2|  1.5|      0.50| 0.000  |
#> | Model 11 |    4|  1.5|      0.50| 0.000  |
#> | Model 12 |    8|  1.5|      0.50| 0.000  |
#> | Model 13 |    1|  2.5|      0.50| 0.000  |
#> | Model 14 |    2|  2.5|      0.50| 0.000  |
#> | Model 15 |    4|  2.5|      0.50| 0.000  |
#> | Model 16 |    8|  2.5|      0.50| 0.401  |
#> | Model 17 |    1|  0.5|      0.75| 0.000  |
#> | Model 18 |    2|  0.5|      0.75| 0.000  |
#> | Model 19 |    4|  0.5|      0.75| 0.000  |
#> | Model 20 |    8|  0.5|      0.75| 0.000  |
#> | Model 21 |    1|  1.0|      0.75| 0.000  |
#> | Model 22 |    2|  1.0|      0.75| 0.000  |
#> | Model 23 |    4|  1.0|      0.75| 0.000  |
#> | Model 24 |    8|  1.0|      0.75| 0.000  |
#> | Model 25 |    1|  1.5|      0.75| 0.000  |
#> | Model 26 |    2|  1.5|      0.75| 0.000  |
#> | Model 27 |    4|  1.5|      0.75| 0.000  |
#> | Model 28 |    8|  1.5|      0.75| 0.000  |
#> | Model 29 |    1|  2.5|      0.75| 0.000  |
#> | Model 30 |    2|  2.5|      0.75| 0.000  |
#> | Model 31 |    4|  2.5|      0.75| 0.471  |
#> | Model 32 |    8|  2.5|      0.75| 0.128  |
#> +----------+-----+-----+----------+--------+
#> 
#> ----------------------------------------
#>  Diagnostics
#> ----------------------------------------
#> Model 31 (stacking weight 0.471):
#>   - smallest relative Cholesky pivot is 8.4e-09: the covariance matrix is
#>     nearly singular and results may lose accuracy (nearly coincident
#>     locations, or a very small decay parameter).
#> 5 other candidate model(s) with stacking weight at most 0.05 have numerical
#> flags; see the 'diagnostics' element.
#> ----------------------------------------
```

**Step 3.** Run
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md)
to obtain posterior samples from the stacked posterior and then analyze
the output.

``` r

post_samps <- stackedSampler(mod1)

post_beta <- post_samps$beta
summary_beta <- t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod1$X.names
print(summary_beta)
#>                   2.5%        50%      97.5%
#> (Intercept) -2.2758080  0.4622453  2.6466488
#> x1          -0.6267972 -0.5155141 -0.4253295
#> x2           0.2114782  0.3083702  0.4063226
```

> **Note:** The following optional steps are only required if interested
> in posterior predictive inference.

**(Optional) Step 4.** Prepare the inputs for posterior predictive
inference.

``` r

sp_pred <- as.matrix(dat_pred[, c("s1", "s2")])
X_new <- cbind(1, dat_pred$x1, dat_pred$x2)
```

**(Optional) Step 5.** Finally, pass the model output through the
function
[`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)

``` r

mod.pred <- posteriorPredict(mod1,
                             coords_new = sp_pred,
                             covars_new = X_new,
                             joint = FALSE)
```

**(Optional) Step 6.** Once samples from the posterior predictive
distribution is obtained, obtain samples from the stacked posterior
predictive distribution using
[`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md).

``` r

postpred_samps <- stackedSampler(mod.pred)
```

Further, we analyze the posterior predictive distribution of the spatial
process against their corresponding true values.

``` r

postpred_z <- postpred_samps$z.pred
post_z_summ <- t(apply(postpred_z, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
z_combn <- data.frame(z = dat_pred$z_true, zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2], zU = post_z_summ[, 3])
ggplot(data = z_combn, aes(x = z)) +
  geom_errorbar(aes(ymin = zL, ymax = zU), color = "skyblue") +
  geom_point(aes(y = zM), color = "darkblue", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "solid") +
  xlab("True z") + ylab("Posterior predictive of z") + theme_bw() +
  theme(panel.background = element_blank(), aspect.ratio = 1)
```

![](spStack_files/figure-html/unnamed-chunk-17-1.png)

## Additional functionalities

We have devised and demonstrated Bayesian predictive stacking to be an
effective tool for estimating spatial/spatial-temporal regression models
and yielding robust predictions for Gaussian as well as non-Gaussian
data. We develop and exploit analytically accessible distribution theory
pertaining to Bayesian analysis of linear mixed model and generalized
linear mixed models that enables us to directly sample from the
posterior distributions. The focus of this package is on effectively
combining inference across different closed-form posterior distributions
by circumventing inference on weakly identified parameters.

The “Technical Overview” vignette provides a comprehensive review of the
conjugate Bayesian hierarchical models used here. The “Spatial Models”
and “Spatial-Temporal Models” explains in detail how to use the
functions for spatial and spatial-temporal regressions, respectively.
Future developments and investigations will consider zero-inflated
non-Gaussian data and adapting to variants of Gaussian process models
that scale inference to massive datasets by circumventing the Cholesky
decomposition of dense covariance matrices.

## References

Pan, Soumyakanti, Lu Zhang, Jonathan R. Bradley, and Sudipto Banerjee.
2025. “Bayesian Inference for Spatial-Temporal Non-Gaussian Data Using
Predictive Stacking.” *Bayesian Analysis* (In press).
<https://doi.org/10.1214/25-BA1582>.

Zhang, Lu, Wenpin Tang, and Sudipto Banerjee. 2025. “Bayesian
Geostatistics Using Predictive Stacking.” *Journal of the American
Statistical Association* (In press) (January).
<https://doi.org/10.1080/01621459.2025.2566449>.
