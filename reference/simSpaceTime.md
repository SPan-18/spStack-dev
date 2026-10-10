# Synthetic point-referenced spatial-temporal data

Dataset of size 500 with spatial coordinates sampled uniformly from the
unit square and temporal coordinates sampled uniformly from the unit
interval, two covariates, an intercept and a slope of `x1` that vary
over space and time, and Gaussian and Poisson responses that share them.

## Usage

``` r
data(simSpaceTime)
```

## Format

a `data.frame` object with 500 rows and columns

- `s1, s2`:

  2-D coordinates in the unit square.

- `t_coords`:

  temporal coordinates in the unit interval.

- `x1, x2`:

  covariates sampled from the standard normal distribution.

- `y_gauss`:

  Gaussian response.

- `y_pois`:

  Poisson count response.

- `z1_true`:

  true spatial-temporal effect associated with the intercept.

- `z2_true`:

  true spatial-temporal effect associated with `x1`.

## Details

With \\\ell = (s, t)\\, the varying intercept is a wave travelling
across space over time and the varying slope of \\x_1\\ changes with
\\s_2\\ and \\t\\, \$\$ z_1(\ell) = \sin\\2 \pi (s_1 - t)\\, \quad
z_2(\ell) = \cos(2 \pi s_2) \cos(\pi t), \$\$ deterministic surfaces
that are not draws from Gaussian processes. The responses are generated
as \$\$ \begin{aligned} y\_{\mathrm{gauss}}(\ell) &\sim N(2 + 5
x_1(\ell) - x_2(\ell) + z_1(\ell) + x_1(\ell) z_2(\ell), 0.5^2),\\
y\_{\mathrm{pois}}(\ell) &\sim \mathrm{Poisson}(\exp\\2 - 0.5
x_1(\ell) + 0.3 x_2(\ell) + z_1(\ell) + x_1(\ell) z_2(\ell)\\).
\end{aligned} \$\$ This data can be generated with the code given in the
example.

## See also

[simSpatial](https://span-18.github.io/spStack-dev/reference/simSpatial.md)

## Author

Soumyakanti Pan <span18@ucla.edu>

## Examples

``` r
set.seed(1726)
n <- 500
s1 <- runif(n)
s2 <- runif(n)
t_coords <- runif(n)
x1 <- rnorm(n)
x2 <- rnorm(n)
z1 <- sin(2 * pi * (s1 - t_coords))
z2 <- cos(2 * pi * s2) * cos(pi * t_coords)
dat <- data.frame(
  s1 = s1, s2 = s2, t_coords = t_coords, x1 = x1, x2 = x2,
  y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z1 + x1 * z2, sd = 0.5),
  y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z1 + x1 * z2)),
  z1_true = z1, z2_true = z2
)
all.equal(dat, simSpaceTime)
#> [1] TRUE
```
