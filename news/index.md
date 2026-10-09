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
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  faster GLM backend. The pre-processing of the conjugate sampler needs
  one n x n Cholesky factorization instead of several, and is about 3x
  faster for
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md).
  Exact leave-one-out and K-fold cross-validation update the full-data
  pre-processing for each held-out site or fold (O(n^2) per site)
  instead of recomputing it (O(n^3)), with a direct recomputation as a
  fallback. Posterior and Monte Carlo draws are projected in blocks of
  64 with level-3 BLAS, drawing the random variates in the original
  order: results are unchanged up to floating-point rounding, and
  speed-ups are largest with an optimized BLAS (e.g. about 3.5x for
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)
  and 5.5x for
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  with OpenBLAS). Results are unchanged up to floating-point rounding.
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  candidate models that differ only in `boundary` share all the
  pre-processing (full data, and each leave-one-out site or
  cross-validation fold): candidates with the same (phi, nu), or (phi_s,
  phi_t), are fitted in one call that loops over their boundary values
  inside each fold. With three boundary values and 10-fold CV this cuts
  the run time by 25-55% with OpenBLAS (e.g. n = 1000: 23 s to 16 s for
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  26 s to 19 s for
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md);
  n = 2000: 146 s to 66 s for
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md))
  and by up to 15% with the reference BLAS. The posterior samples of
  these candidates are now drawn before their leave-one-out draws, so
  results for a given seed differ from earlier versions at the Monte
  Carlo level; a candidate fitted on its own (or a group with a single
  boundary value) is unchanged.
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  new advanced tag `CV.update` in `loopd.controls` (`'auto'`, `'update'`
  or `'direct'`) choosing how the pre-processing of each
  cross-validation fold is obtained: by deletion updates of the
  full-data factors (scalar loops; faster with the reference BLAS) or by
  recomputing it on the fold (level-3 BLAS; faster with an optimized
  BLAS, e.g. a further 25-50% for
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
  with OpenBLAS at n = 1000-2000). The default `'auto'` uses `'direct'`
  when R reports an optimized BLAS (OpenBLAS, MKL, BLIS,
  Accelerate/vecLib, ATLAS). Results agree up to floating-point
  rounding.
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)
  and
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  use `'auto'`.
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  leaving out `loopd.controls`, or the `CV.K` or `nMC` tag
  ([`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)),
  stopped with an error; the documented defaults are now used.
- GLMs: the latent pseudo-data are now drawn on the log scale, as log
  G1 - log G2 with gamma variates G1, G2 (binomial, binary) and as log G
  with an underflow-safe draw for shape \< 1 (Poisson). Previously, a
  beta draw that rounded to 1, or a gamma draw that underflowed to 0,
  produced infinite or NaN samples and leave-one-out predictive
  densities; this happened with positive probability for small
  `boundary` (e.g. about 1 in 2 million draws at `boundary = 0.4` for
  binary data). The sampled distribution is unchanged, but results for a
  given seed differ from earlier versions at the Monte Carlo level
  (binomial and binary data; Poisson data with zero counts).
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md):
  `boundary` is no longer raised to 0.4 for binary data. A message is
  given when `boundary < 0.1`.
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  each candidate `boundary` must lie in (0, 1), and a message is given
  when any is below 0.1.
- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  fixed a heap buffer overflow in K-fold cross-validation when the
  number of varying coefficients exceeded the number of fixed-effect
  covariates.
- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  `process.type` is now checked to be one of ‘independent’,
  ‘independent.shared’ or ‘multivariate’. The undocumented, incomplete
  “multivariate2” option has been removed.
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  and
  [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md)
  for GLMs: faster and leaner (the kriging means of 64 draws at a time
  with level-3 BLAS, one n x n matrix fewer, the joint conditional
  covariance built only for joint prediction). The C++ code now saves
  and restores R’s random number state, so `.Random.seed` advances and
  the draws no longer overlap with later random numbers in the session;
  pointwise prediction variances are clamped at 0 against rounding.
- [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md)
  for
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)/[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)
  fits with `process.type = "independent"`: the posterior of each
  process’s scale used the shape (nu.z + n\*r)/2 of the shared-scale
  model; it is now (nu.z + n)/2. Recovered scales from earlier versions
  were too concentrated and biased downwards.
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  for
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)/[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)
  fits with `process.type = "multivariate"` and `joint = FALSE`:
  `mu.pred` and `y.pred` were computed from the unscaled noise instead
  of the predicted `z.pred`; fixed.
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
