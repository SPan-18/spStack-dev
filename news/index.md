# Changelog

## spStack (development version)

- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
  [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md):
  faster, leaner Gaussian backend. Exact leave-one-out predictive
  densities are computed in closed form, the fit holds a single n x n
  matrix, distances are computed from `coords` in C++, Matern
  correlations use closed forms for nu = 0.5, 1.5, 2.5, and
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)
  builds each correlation matrix once per (phi, nu). Results are
  unchanged up to floating-point rounding.
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md):
  distances are computed from `coords` in C++, and
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
  builds each correlation matrix once per (phi, nu).
- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md):
  the posterior scale of the variance is computed in residual form,
  avoiding loss of precision when the signal is large relative to the
  noise; Cholesky failures now stop with an informative error.
- [`cholUpdateDel()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md),
  [`cholUpdateDelBlock()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md):
  indices are now checked to be single integers between 1 and `n`; an
  index of 0 previously returned a zero matrix or crashed R.
- All model-fitting functions now stop with an informative error if
  coordinates are duplicated; for
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  and
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md),
  a duplicate must coincide in both space and time.
- PSIS leave-one-out predictive densities (`loopd.method = "PSIS"`) are
  reimplemented in C++ following the loo package, with O(S) memory;
  Pareto k diagnostics are returned as `loopd.pareto_k`.
- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md):
  the inverse-gamma prior is now placed on the measurement error
  variance `sigmaSq`; posterior samples of the spatial variance are
  returned as `sigmaSq.z`. The default prior is now `priors = "flat"`,
  i.e., p(beta, sigmaSq) proportional to 1/sigmaSq.

## spStack 1.1.3

CRAN release: 2026-03-08

- Documentation Update: Pan, Zhang, Bradley and Banerjee (2025) accepted
  at Bayesian Analysis.
- Migrate to CVXR 1.8.1 API and resolve CRAN Results errors.
- Add a fallback option to the optimization routine.

## spStack 1.1.2

CRAN release: 2025-10-04

- Documentation Update: Zhang, Tang and Banerjee (2025) accepted at
  JASA.

## spStack 1.1.1

CRAN release: 2025-07-14

- `lmulm_XTilde_VC()`, `lmulv_XTilde_VC()`: Fixed address sanitizer
  issue with string comparison with pointer to string literal.

## spStack 1.1.0

CRAN release: 2025-07-12

- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  New functions for spatially-temporally varying coefficients GLM.

- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md):
  New functions for posterior predictive inference using predictive
  stacking.

- [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md):
  Utility for recovering posterior samples of scale parameters in
  spatial and spatial-temporal GLMs.

## spStack 1.0.1

CRAN release: 2024-10-08

- Fix a memory leak issue.

## spStack 1.0.0

CRAN release: 2024-10-03

- Initial CRAN submission.
