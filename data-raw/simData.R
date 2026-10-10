## Code to generate the synthetic datasets shipped with spStack: simSpatial and
## simSpaceTime. Run from the package root: source("data-raw/simData.R").

## simSpatial: point-referenced data on the unit square sharing the covariates
## x1, x2 and the spatial effect z, with one response for each family.
set.seed(1729)
n <- 500
s1 <- runif(n)
s2 <- runif(n)
x1 <- rnorm(n)
x2 <- rnorm(n)
z <- 1.5 * sin(4 * pi * sqrt((s1 - 0.5)^2 + (s2 - 0.5)^2))     # ripple
n_trials <- sample(5:20, n, replace = TRUE)
simSpatial <- data.frame(
  s1 = s1, s2 = s2, x1 = x1, x2 = x2,
  y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z, sd = 0.5),
  y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z)),
  y_binom = rbinom(n, n_trials, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  n_trials = n_trials,
  y_bin = rbinom(n, 1, plogis(0.5 - 0.5 * x1 + 0.5 * x2 + z)),
  y_svc = rnorm(n, 2 + (5 + z) * x1 - x2, sd = 0.5),
  z_true = z
)

## simSpaceTime: spatial-temporal data with a varying intercept (a wave
## travelling across space over time) and a varying slope of x1.
set.seed(1726)
n <- 500
s1 <- runif(n)
s2 <- runif(n)
t_coords <- runif(n)
x1 <- rnorm(n)
x2 <- rnorm(n)
z1 <- sin(2 * pi * (s1 - t_coords))                              # varying intercept
z2 <- cos(2 * pi * s2) * cos(pi * t_coords)                # varying slope of x1
simSpaceTime <- data.frame(
  s1 = s1, s2 = s2, t_coords = t_coords, x1 = x1, x2 = x2,
  y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z1 + x1 * z2, sd = 0.5),
  y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z1 + x1 * z2)),
  z1_true = z1, z2_true = z2
)

save(simSpatial, file = "data/simSpatial.rda", compress = "xz")
save(simSpaceTime, file = "data/simSpaceTime.rda", compress = "xz")
