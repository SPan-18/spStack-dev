# Synthetic point-referenced spatial-temporal Poisson count data simulated using spatially-temporally varying coefficients

Dataset of size 500, with a Poisson distributed response variable
indexed by spatial and temporal coordinates sampled uniformly from the
unit square. The model includes an intercept and a covariate with
spatially and temporally varying coefficients, and spatial-temporal
random effects induced by a Matérn covariogram.

## Usage

``` r
data(sim_stvcPoisson)
```

## Format

a `data.frame` object.

- `s1, s2`:

  2-D coordinates; latitude and longitude.

- `t_coords`:

  temporal coordinates.

- `x1`:

  a covariate sampled from the standard normal distribution.

- `y`:

  response vector.

- `z1_true`:

  true spatial-temporal random effect associated with the intercept.

- `z2_true`:

  true spatial-temporal random effect associated with the covariate.

## Details

With \\n = 500\\, the count data is simulated using \$\$ \begin{aligned}
y(s_i) &\sim \mathrm{Poisson}(\lambda(s_i)), i = 1, \ldots, n,\\ \log
\lambda(s_i) &= x(s_i)^\top \beta + x(s_i)^\top z(s_i), \end{aligned}
\$\$ where the spatial-temporal random effects \\z(s) = (z_1(s),
z_2(s))^\top\\ with independent processes \\z_j(s) \sim GP(0, \sigma_j^2
R(\cdot, \cdot; \theta_j))\\ for \\j = 1, 2\\, and with \\R(\cdot,
\cdot; \theta_j)\\ given by \$\$ R((s, t), (s', t'); \theta_j) =
\frac{1}{\phi\_{tj} \|t - t'\|^2 + 1} \exp \left( - \frac{\phi\_{sj}
\lVert s - s' \rVert}{\sqrt{1 + \phi\_{tj} \|t - t'\|^2}} \right) , \$\$
where \\\phi_s\\ is the spatial decay parameter, \\\phi_t\\ is the
temporal decay parameter. We have sampled the data with \\\beta = (2,
-0.5)\\, \\\phi\_{s1} = 2\\, \\\phi\_{s2} = 3\\, \\\phi\_{t1} = 2\\,
\\\phi\_{t2} = 4\\, \\\sigma^2_1 = \sigma^2_2 = 1\\. This data can be
generated with the code as given in the example below.

## See also

[simPoisson](https://span-18.github.io/spStack-dev/reference/simPoisson.md),
[simGaussian](https://span-18.github.io/spStack-dev/reference/simGaussian.md),
[simBinom](https://span-18.github.io/spStack-dev/reference/simBinom.md),
[simBinary](https://span-18.github.io/spStack-dev/reference/simBinary.md)

## Examples

``` r
rmvn <- function(n, mu = 0, V = matrix(1)) {
p <- length(mu)
  if (any(is.na(match(dim(V), p))))
    stop("error: dimension mismatch.")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}
set.seed(1726)
n <- 500
beta <- c(2, -0.5)
p <- length(beta)
X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
X_tilde <- X
phi_s <- c(2, 3)
phi_t <- c(2, 4)
S <- data.frame(s1 = runif(n, 0, 1), s2 = runif(n, 0, 1))
Tm <- runif(n)
dist_S <- as.matrix(dist(as.matrix(S)))
dist_T <- as.matrix(dist(as.matrix(Tm)))
Vz1 <- 1/(1 + phi_t[1] * dist_T^2) * exp(- (phi_s[1] * dist_S) / sqrt(1 + phi_t[1] * dist_T^2))
Vz2 <- 1/(1 + phi_t[2] * dist_T^2) * exp(- (phi_s[2] * dist_S) / sqrt(1 + phi_t[2] * dist_T^2))
z1 <- rmvn(1, rep(0, n), Vz1)
z2 <- rmvn(1, rep(0, n), Vz2)
muFixed <- X %*% beta
muSpT <- X_tilde[, 1] * z1 + X_tilde[, 2] * z2
mu <- muFixed + muSpT
y <- sapply(1:n, function(x) rpois(n = 1, lambda = exp(mu[x])))
dat <- cbind(S, Tm, X[, -1], y, z1, z2)
names(dat) <- c("s1", "s2", "t_coords", paste("x", 1:(p - 1), sep = ""), "y", "z1_true", "z2_true")
```
