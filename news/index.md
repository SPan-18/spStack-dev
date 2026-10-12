# Changelog

## spStack (development version)

- New model families:
  [`svcLMexact()`](https://span-18.github.io/spStack-dev/reference/svcLMexact.md),
  [`svcLMstack()`](https://span-18.github.io/spStack-dev/reference/svcLMstack.md),
  [`svcGLMexact()`](https://span-18.github.io/spStack-dev/reference/svcGLMexact.md)
  and
  [`svcGLMstack()`](https://span-18.github.io/spStack-dev/reference/svcGLMstack.md)
  for spatially varying coefficients model, and
  [`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md),
  [`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md)
  for Gaussian spatially-temporally varying coefficients model.
- Faster, leaner and more robust backend: closed-form exact LOO for the
  Gaussian models, deletion updates for GLM LOO and K-fold CV, PSIS in
  C++, informative errors, and a new `diagnostics` element.
- Breaking changes: in
  [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)/[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)
  the inverse-gamma prior is on the noise variance `sigmaSq`; stacked
  fits store `model.params` and `stacking.summary`.
- New lean synthetic datasets `simSpatial` and `simSpaceTime` replace
  the old datasets and `sim_spData()`; all examples and vignettes use
  them.
  [`surfaceplot()`](https://span-18.github.io/spStack-dev/reference/surfaceplot.md)/[`surfaceplot2()`](https://span-18.github.io/spStack-dev/reference/surfaceplot2.md)
  use a colorblind-friendly palette.

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
