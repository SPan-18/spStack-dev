# Package index

## Stacking

Implements Bayesian predictive stacking for spatial/spatial-temporal
linear and generalized linear models.

- [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)
  : Bayesian spatial linear model using predictive stacking
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
  : Bayesian spatial generalized linear model using predictive stacking
- [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)
  : Bayesian spatially-temporally varying coefficients generalized
  linear model using predictive stacking
- [`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md)
  : Bayesian spatially-temporally varying coefficients linear model
  using predictive stacking
- [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md)
  : Create a collection of candidate models for stacking

## Exact samplers and finds LOO-PD

Samples from the posterior distribution conditional on fixed values of
spatial/spatial-temporal process parameters and finds leave-one-out
predictive densities

- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)
  : Univariate Bayesian spatial linear model
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)
  : Univariate Bayesian spatial generalized linear model
- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  : Bayesian spatially-temporally varying generalized linear model
- [`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md)
  : Bayesian spatially-temporally varying coefficients linear model

## Downstream posterior sampling

Posterior recovery of scale parameters and posterior predictive
inference

- [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md)
  : Recover posterior samples of scale parameters of
  spatial/spatial-temporal generalized linear models
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  : Prediction of latent process at new spatial or temporal locations

## Analysis of output

### Stacked posterior

Helps analyze output of stacking algorithms

- [`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md)
  : Sample from the stacked posterior distribution
- [`get_stacking_weights()`](https://span-18.github.io/spStack-dev/reference/get_stacking_weights.md)
  : Optimal stacking weights

### Surface plots

Creates and plots interpolated spatial surfaces

- [`surfaceplot()`](https://span-18.github.io/spStack-dev/reference/surfaceplot.md)
  : Make a surface plot
- [`surfaceplot2()`](https://span-18.github.io/spStack-dev/reference/surfaceplot2.md)
  : Make two side-by-side surface plots

## Matrix algorithms

Provides compiled code for various Cholesky factor updates helpful for
those designing their own cross-validation algorithms

- [`cholUpdateRankOne()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md)
  [`cholUpdateDel()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md)
  [`cholUpdateDelBlock()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md)
  : Different Cholesky factor updates

## Synthetic datasets

Simulated spatial and spatial-temporal data on the unit square

- [`simSpatial`](https://span-18.github.io/spStack-dev/reference/simSpatial.md)
  : Synthetic point-referenced spatial data
- [`simSpaceTime`](https://span-18.github.io/spStack-dev/reference/simSpaceTime.md)
  : Synthetic point-referenced spatial-temporal data

## Other utilities

- [`spStack`](https://span-18.github.io/spStack-dev/reference/spStack-package.md)
  [`spStack-package`](https://span-18.github.io/spStack-dev/reference/spStack-package.md)
  : spStack: Bayesian Geostatistics Using Predictive Stacking
- [`iDist()`](https://span-18.github.io/spStack-dev/reference/iDist.md)
  : Calculate distance matrix
