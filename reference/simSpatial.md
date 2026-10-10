# Synthetic point-referenced spatial data

Dataset of size 500 with spatial coordinates sampled uniformly from the
unit square, two covariates, a spatial effect with a ripple pattern and
one response for each of the Gaussian, Poisson, binomial and binary
families, plus a Gaussian response with a spatially varying slope. All
the responses share the covariates and the spatial effect, so that every
spatial model in the package can be fitted to the same data.

## Usage

``` r
data(simSpatial)
```

## Format

a `data.frame` object with 500 rows and columns

- `s1, s2`:

  2-D coordinates in the unit square.

- `x1, x2`:

  covariates sampled from the standard normal distribution.

- `y_gauss`:

  Gaussian response.

- `y_pois`:

  Poisson count response.

- `y_binom`:

  binomial count response, out of `n_trials` trials.

- `n_trials`:

  number of trials of the binomial response (5 to 20).

- `y_bin`:

  binary response.

- `y_svc`:

  Gaussian response in which `z` is the spatially varying part of the
  slope of `x1`.

- `z_true`:

  true spatial effect that generated the data.

## Details

The spatial effect is the ripple \$\$ z(s) = 1.5 \sin(4 \pi \lVert s -
(0.5, 0.5) \rVert), \$\$ a deterministic surface that is not a draw from
a Gaussian process, and the responses are generated as \$\$
\begin{aligned} y\_{\mathrm{gauss}}(s) &\sim N(2 + 5 x_1(s) - x_2(s) +
z(s), 0.5^2),\\ y\_{\mathrm{pois}}(s) &\sim \mathrm{Poisson}(\exp\\2 -
0.5 x_1(s) + 0.3 x_2(s) + z(s)\\),\\ y\_{\mathrm{binom}}(s) &\sim
\mathrm{Binomial}(m(s), \mathrm{ilogit}\\0.5 - 0.5 x_1(s) + 0.5 x_2(s) +
z(s)\\),\\ y\_{\mathrm{bin}}(s) &\sim
\mathrm{Bernoulli}(\mathrm{ilogit}\\0.5 - 0.5 x_1(s) + 0.5 x_2(s) +
z(s)\\),\\ y\_{\mathrm{svc}}(s) &\sim N(2 + \\5 + z(s)\\ x_1(s) -
x_2(s), 0.5^2), \end{aligned} \$\$ where \\m(s)\\ is the number of
trials. This data can be generated with the code given in the example.

## See also

[simSpaceTime](https://span-18.github.io/spStack-dev/reference/simSpaceTime.md)

## Author

Soumyakanti Pan <span18@ucla.edu>

## Examples

``` r
set.seed(1729)
n <- 500
s1 <- runif(n)
s2 <- runif(n)
x1 <- rnorm(n)
x2 <- rnorm(n)
z <- 1.5 * sin(4 * pi * sqrt((s1 - 0.5)^2 + (s2 - 0.5)^2))
n_trials <- sample(5:20, n, replace = TRUE)
dat <- data.frame(
  s1 = s1, s2 = s2, x1 = x1, x2 = x2,
  y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z, sd = 0.5),
  y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z)),
  y_binom = rbinom(n, n_trials, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  n_trials = n_trials,
  y_bin = rbinom(n, 1, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  y_svc = rnorm(n, 2 + (5 + z) * x1 - x2, sd = 0.5),
  z_true = z
)
all.equal(dat, simSpatial)
#> [1] TRUE
```
