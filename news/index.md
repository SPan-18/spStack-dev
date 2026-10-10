# Changelog

## spStack (development version)

- New synthetic datasets `simSpatial` and `simSpaceTime` replace
  `simGaussian`, `simPoisson`, `simBinom`, `simBinary` and
  `sim_stvcPoisson`, and `sim_spData()` is removed. `simSpatial` has one
  set of locations, covariates `x1`, `x2` and a spatial effect with a
  ripple pattern, shared by Gaussian, Poisson, binomial and binary
  responses (and a Gaussian response with a spatially varying slope);
  `simSpaceTime` has a varying intercept that travels across space over
  time and a varying slope of `x1`, shared by Gaussian and Poisson
  responses. The true surfaces are deterministic, so the examples and
  vignettes show how well the Gaussian process models recover them. All
  examples and vignettes now use these two datasets, and the examples no
  longer draw plots. The code that generates the data is in `data-raw/`
  and in the examples of the dataset help pages. Code that loads the old
  datasets must be updated.
- [`surfaceplot()`](https://span-18.github.io/spStack-dev/reference/surfaceplot.md),
  [`surfaceplot2()`](https://span-18.github.io/spStack-dev/reference/surfaceplot2.md):
  the default palette is now the colorblind-friendly diverging
  ColorBrewer palette ‘RdBu’ (it was ‘RdYlBu’).
- [`stvcLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcLMexact.md),
  [`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md):
  new functions for the Bayesian linear model with spatially-temporally
  varying coefficients (Gaussian response). Each of the r varying
  coefficients has an independent spatial-temporal process (Gneiting
  correlation) with its own (`process.type = "independent"`) or a common
  (`"independent.shared"`) set of decay parameters `phi_s`, `phi_t` and
  noise-to-spatial variance ratio `noise_sp_ratio`; the inverse-gamma
  prior is on the noise variance `sigmaSq`, and the process variances
  are returned as `sigmaSq.z = sigmaSq / noise_sp_ratio`. The joint
  posterior is sampled exactly by composition; the nr-dimensional latent
  process is drawn with the update of Bhattacharya, Chakraborty and
  Mallick (2016) (Matheron’s rule), which needs only n x n Cholesky
  factorizations and processes the draws with level-3 BLAS (1.3-5.9x
  faster than factorizing its nr x nr posterior covariance for r \>= 2,
  n \>= 1000). Exact leave-one-out predictive densities are computed in
  closed form from the same factorizations (3-7% added to a fit), or by
  PSIS.
  [`stvcLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcLMstack.md)
  stacks candidate models given by
  [`candidateModels()`](https://span-18.github.io/spStack-dev/reference/candidateModels.md),
  building and factorizing each correlation matrix once per distinct
  (`phi_s`, `phi_t`).
  [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  and
  [`stackedSampler()`](https://span-18.github.io/spStack-dev/reference/stackedSampler.md)
  support both.
- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md):
  with `loopd = TRUE` and `loopd.method` not supplied, the call stopped
  with “loopd.method must be specified” although the argument has a
  default (`"exact"`); the default is now used.
- [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md):
  without `loopd.method`, the message said ‘exact’ would be used but the
  argument was never set, so the call stopped; it now uses ‘exact’.
- [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  candidate models for `process.type = "independent"` were required to
  have `phi_s` and `phi_t` of length 2 whatever the number of varying
  coefficients; they must now have length r.
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  for the spatial-temporal models: a `coords_new` or `covars_new` list
  without the required tags (`sp` and `time`, `fixed` and `vc`) is now
  reported with an informative error.
- Stacking functions with `process.type = "independent"` and PSIS: the
  Pareto k message of a candidate model in the “Diagnostics” section was
  repeated once per process; it is now printed once.
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
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md),
  [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md):
  a failed Cholesky factorization (e.g. a numerically singular
  correlation matrix from nearly coincident locations or very strong
  correlation, a prior covariance `V.beta` or `iw.scale` that is not
  positive definite) now stops with an informative error. Previously a
  message was printed to the console and the computation continued with
  invalid factors. In leave-one-out and cross-validation, a failed
  deletion update still falls back to the direct computation, and only a
  failure of that stops.
- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  with `loopd = TRUE`, and
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md),
  for `process.type = "multivariate"`: the prior draw of the latent
  process in the posterior sampler used the Cholesky factor of the
  inverse-Wishart draw without zeroing its upper triangle, so the
  returned posterior samples of `beta`, `z` and `xi` had the wrong
  distribution. Fixed; the samples now match
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
  with `loopd = FALSE` for the same seed. Leave-one-out predictive
  densities and stacking weights were not affected.
- [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  K-fold cross-validation now saves R’s random number state at the end,
  so later random numbers in the session no longer repeat the ones used
  by the cross-validation draws.
- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  new `diagnostics` element collecting the fit diagnostics, with no
  warnings issued (for the GLMs the pivot is that of chol(Vz), and for
  `process.type = "independent"` there is one row per process).
  `diagnostics$numerical` (one row per model) holds the smallest
  relative Cholesky pivot of the n x n factorizations and the
  correlations of the two farthest-apart and the two closest locations,
  taken from quantities the fit computes anyway (no extra
  factorization). `diagnostics$pareto` holds the Pareto k values of PSIS
  leave-one-out (previously `loopd.pareto_k` and a warning), with the
  threshold and the number of values above it. The stacking functions:
  `diagnostics$solver` holds the solver details of the stacking weights
  (solver used and status, installed and requested solvers, search
  order, each attempt with its status, and whether the
  [`loo::stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  fallback was used); this replaces the top-level `solver` and
  `solver.status` elements and the solver messages printed with
  `verbose = TRUE`. With `verbose = TRUE`, a “Diagnostics” section is
  printed only if there is an issue: a nearly singular covariance matrix
  (pivot \< 1e-8), an effective range far beyond the extent of the data
  (farthest-pair correlation \> 0.95), nearly uncorrelated locations
  (closest-pair correlation \< 0.05), Pareto k values above the
  threshold, or a solver problem (requested solver not installed,
  inaccurate solution, fallback). In the stacking functions, numerical
  flags are detailed only for candidates with stacking weight above
  0.05.
- [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  new element `model.params`, a list with the parameters of each
  candidate model as a named list (`phi`, `nu`, `noise_sp_ratio` or
  `boundary`; `phi_s`, `phi_t`, `boundary` for
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)).
  [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  and
  [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md)
  now read the model parameters only from it (by name, through one
  internal accessor), instead of from the printed table of each model
  class. The table of candidate models and stacking weights is now
  `stacking.summary` in all three functions, and the `candidate.models`
  element is removed (it was this table in
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)/[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
  and the parameter list in
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)).
  Stacked fits saved with earlier versions must be refitted to be used
  with
  [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md)
  or
  [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.md).
- [`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md):
  for the exponential correlation function `nu` is now `NA` (it was 0)
  in `model.params`.
- [`get_stacking_weights()`](https://span-18.github.io/spStack-dev/reference/get_stacking_weights.md):
  returns the solver details in a new `details` element; the error
  message when both CVXR and the loo fallback fail now reports the last
  CVXR error (it was always NULL).
- [`get_stacking_weights()`](https://span-18.github.io/spStack-dev/reference/get_stacking_weights.md):
  the
  [`loo::stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  fallback (used when no CVXR solver reaches an optimal solution) was
  given the shifted predictive densities instead of the log predictive
  densities, so its weights were not the stacking optimum; it now gets
  the log densities.
- `loo` moves from Imports to Suggests: it is used only by that
  fallback. If it is needed but not installed, the stacking functions
  return the fitted models with `NA` stacking weights and a message
  giving the code that computes the weights once `loo` is installed.
- All model-fitting functions: missing values in the response,
  covariates, binomial trials or coordinates now stop with an
  informative error. Previously
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) silently
  dropped incomplete rows of the data but not of the coordinates, which
  stopped with a misleading error about the number of coordinate rows.
- All model-fitting functions: `beta.norm[[2]]`
  ([`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md),
  [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)),
  `V.beta` and `IW.scale` (GLMs) are checked to be symmetric and
  positive definite; a non-symmetric matrix was previously used through
  its lower triangle without notice.
- All model-fitting functions: arguments passed through `...` are
  reported with a warning (they are not used). The previous check
  compared them with the arguments of the calling function, so it could
  warn or stay silent incorrectly when the function was called from
  another function.
- [`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md),
  [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md),
  [`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md),
  [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md):
  exact leave-one-out and K-fold cross-validation can be interrupted by
  the user (Esc / Ctrl-C) between sites or folds, with all memory
  released.
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.md):
  the check that new locations differ from the observed ones bins the
  coordinates instead of comparing all pairs (5000 new and 5000 observed
  locations: about 38 s to 0.1 s), with the same result.
- [`cholUpdateRankOne()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md),
  [`cholUpdateDel()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md),
  [`cholUpdateDelBlock()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.md)
  with `lower = FALSE`: the C++ code now works on a copy, so the input
  matrix can never be modified in place.
- Internal: unused C++ code removed (the superseded one-draw-at-a-time
  GLM projections, the unreachable GLM PSIS branch and its helpers, and
  other dead routines); the unused `t(X) X` and `t(XTilde) X` products
  in the varying-coefficients pre-processing are no longer computed.
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
  Pareto k diagnostics are returned in the `diagnostics` element (see
  below).
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
