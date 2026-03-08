# Simulate spatial data on unit square

Generates synthetic spatial data of different types where the spatial
co-ordinates are sampled uniformly on an unit square. Different types
include point-referenced Gaussian, Poisson, binomial and binary data.
The design includes an intercept and fixed covariates sampled from a
standard normal distribution.

## Usage

``` r
sim_spData(n, beta, cor.fn, spParams, spvar, deltasq, family, n_binom)
```

## Arguments

- n:

  sample size.

- beta:

  a \\p\\-dimensional vector of fixed effects.

- cor.fn:

  a quoted keyword that specifies the correlation function used to model
  the spatial dependence structure among the observations. Supported
  covariance model key words are: `'exponential'` and `'matern'`.

- spParams:

  a numeric vector containing spatial process parameters - e.g., spatial
  decay and smoothness.

- spvar:

  value of spatial variance parameter.

- deltasq:

  value of noise-to-spatial variance ratio.

- family:

  a character specifying the distribution of the response as a member of
  the exponential family. Valid inputs are `'gaussian'`, `'poisson'`,
  `'binary'`, and `'binomial'`.

- n_binom:

  necessary only when `family = 'binomial'`. Must be a vector of length
  `n` that will specify the number of trials for each observation. If it
  is of length 1, then that value is considered to be the common value
  for the number of trials for all `n` observations.

## Value

a `data.frame` object containing the columns -

- `s1, s2`:

  2D-coordinates in unit square

- `x1, x2, ...`:

  covariates, not including intercept

- `y`:

  response

- `n_trials`:

  present only when binomial data is generated

- `z_true`:

  true spatial effects with which the data is generated

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1729)
n <- 10
beta <- c(2, 5)
phi0 <- 2
nu0 <- 0.5
spParams <- c(phi0, nu0)
spvar <- 0.4
deltasq <- 1
sim1 <- sim_spData(n = n, beta = beta, cor.fn = "matern",
                   spParams = spParams, spvar = spvar, deltasq = deltasq,
                   family = "gaussian")
```
