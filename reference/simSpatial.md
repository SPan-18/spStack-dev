# Synthetic point-referenced spatial data

Dataset of size 500 with spatial coordinates sampled uniformly from the
unit square, two covariates, a spatial effect with a checkerboard
pattern and one response for each of the Gaussian, Poisson, binomial and
binary families, plus a Gaussian and a Poisson response with spatially
varying coefficients (a varying intercept with a ring pattern and a
varying slope of `x1`). All the responses share the covariates and the
spatial effects, so that every spatial model in the package can be
fitted to the same data.

## Usage

``` r
data(simSpatial)
```

## Format

a `data.frame` object with 500 rows and 13 columns

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

  Gaussian response with the spatially varying intercept `z2` and the
  spatially varying slope `5 + z` of `x1`.

- `y_svc_pois`:

  Poisson count response with the spatially varying intercept `z2` and
  the spatially varying slope `-0.5 + z` of `x1`.

- `z_true`:

  true spatial effect (checkerboard) that generated the data; the
  spatially varying part of the slope of `x1` in `y_svc` and
  `y_svc_pois`.

- `z2_true`:

  true spatially varying intercept (rings) of `y_svc` and `y_svc_pois`.

## Details

The spatial effect is the checkerboard \$\$ z(s) = 1.5 \sin(2 \pi s_1)
\sin(2 \pi s_2), \$\$ and the spatially varying intercept is the rings
\$\$ z_2(s) = \cos(4 \pi \lVert s - s_0 \rVert) - c, \quad s_0 = (0.5,
0.5), \$\$ where the constant \\c \approx 0.129\\ makes the mean of
\\z_2\\ over the unit square zero (the mean of \\z\\ is zero by
symmetry, and the two surfaces are orthogonal on the unit square). Both
are deterministic surfaces, not draws from a Gaussian process. The
responses are generated as \$\$ \begin{aligned} y\_{\mathrm{gauss}}(s)
&\sim N(2 + 5 x_1(s) - x_2(s) + z(s), 0.5^2),\\ y\_{\mathrm{pois}}(s)
&\sim \mathrm{Poisson}(\exp\\2 - 0.5 x_1(s) + 0.3 x_2(s) + z(s)\\),\\
y\_{\mathrm{binom}}(s) &\sim \mathrm{Binomial}(m(s),
\mathrm{ilogit}\\0.5 - 0.5 x_1(s) + 0.5 x_2(s) + z(s)\\),\\
y\_{\mathrm{bin}}(s) &\sim \mathrm{Bernoulli}(\mathrm{ilogit}\\0.5 - 0.5
x_1(s) + 0.5 x_2(s) + z(s)\\),\\ y\_{\mathrm{svc}}(s) &\sim N(2 +
z_2(s) + \\5 + z(s)\\ x_1(s) - x_2(s), 0.5^2),\\
y\_{\mathrm{svc\\pois}}(s) &\sim \mathrm{Poisson}(\exp\\2 + z_2(s) +
\\-0.5 + z(s)\\ x_1(s) + 0.3 x_2(s)\\), \end{aligned} \$\$ where
\\m(s)\\ is the number of trials. This data can be generated with the
code given in the example.

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
z <- 1.5 * sin(2 * pi * s1) * sin(2 * pi * s2)
ring <- function(s1, s2) cos(4 * pi * sqrt((s1 - 0.5)^2 + (s2 - 0.5)^2))
g <- (seq_len(1000) - 0.5) / 1000
z2 <- ring(s1, s2) - mean(outer(g, g, ring))
n_trials <- sample(5:20, n, replace = TRUE)
dat <- data.frame(
  s1 = s1, s2 = s2, x1 = x1, x2 = x2,
  y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z, sd = 0.5),
  y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z)),
  y_binom = rbinom(n, n_trials, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  n_trials = n_trials,
  y_bin = rbinom(n, 1, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  y_svc = rnorm(n, 2 + z2 + (5 + z) * x1 - x2, sd = 0.5),
  y_svc_pois = rpois(n, exp(2 + z2 + (-0.5 + z) * x1 + 0.3 * x2)),
  z_true = z,
  z2_true = z2
)
all.equal(dat, simSpatial)
#> [1] TRUE
```
