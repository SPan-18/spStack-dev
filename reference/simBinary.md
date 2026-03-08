# Synthetic point-referenced binary data

Dataset of size 500, with a binary response variable indexed by spatial
coordinates sampled uniformly from the unit square. The model includes
one covariate and spatial random effects induced by a Matérn
covariogram.

## Usage

``` r
data(simBinary)
```

## Format

a `data.frame` object.

- `s1, s2`:

  2-D coordinates; latitude and longitude.

- `x1`:

  a covariate sampled from the standard normal distribution.

- `y`:

  response vector (0/1).

- `z_true`:

  true spatial random effects that generated the data.

## Details

With \\n = 500\\, the binary data is simulated using \$\$
\begin{aligned} y(s_i) &\sim \mathrm{Bernoulli}(\pi(s_i)), i = 1,
\ldots, n,\\ \pi(s_i) &= \mathrm{ilogit}(x(s_i)^\top \beta + z(s_i))
\end{aligned} \$\$ where the function \\\mathrm{ilogit}\\ refers to the
inverse-logit function, the spatial effects \\z \sim N(0, \sigma^2 R)\\
with \\R\\ being a \\n \times n\\ correlation matrix given by the Matérn
covariogram \$\$ R(s, s') = \frac{(\phi \|s-s'\|)^\nu}{\Gamma(\nu)
2^{\nu - 1}} K\_\nu(\phi \|s-s'\|), \$\$ where \\\phi\\ is the spatial
decay parameter and \\\nu\\ the spatial smoothness parameter. We have
sampled the data with \\\beta = (0.5, -0.5)\\, \\\phi = 5\\, \\\nu =
0.5\\, and \\\sigma^2 = 0.4\\. This data can be generated with the code
as given in the example below.

## See also

[simGaussian](https://span-18.github.io/spStack-dev/reference/simGaussian.md),
[simPoisson](https://span-18.github.io/spStack-dev/reference/simPoisson.md),
[simBinom](https://span-18.github.io/spStack-dev/reference/simBinom.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
set.seed(1729)
n <- 500
beta <- c(0.5, -0.5)
phi0 <- 5
nu0 <- 0.5
spParams <- c(phi0, nu0)
spvar <- 0.4
sim1 <- sim_spData(n = n, beta = beta, cor.fn = "matern",
                   spParams = spParams, spvar = spvar, deltasq = deltasq,
                   family = "binary")
plot1 <- surfaceplot(sim1, coords_name = c("s1", "s2"), var_name = "z_true")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the spStack package.
#>   Please report the issue at <https://github.com/SPan-18/spStack-dev/issues>.

library(ggplot2)
plot2 <- ggplot(sim1, aes(x = s1, y = s2)) +
  geom_point(aes(color = factor(y)), alpha = 0.75) +
  scale_color_manual(values = c("red", "blue"), labels = c("0", "1")) +
  guides(alpha = 'none') +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 0.25),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, hjust = 0.25),
        legend.box.just = "center", aspect.ratio = 1)
```
