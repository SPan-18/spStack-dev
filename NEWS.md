# spStack (development version)

* `spLMexact()`, `spLMstack()`, `posteriorPredict()`: faster, leaner Gaussian backend. Exact leave-one-out predictive densities are computed in closed form, the fit holds a single n x n matrix, distances are computed from `coords` in C++, Matern correlations use closed forms for nu = 0.5, 1.5, 2.5, and `spLMstack()` builds each correlation matrix once per (phi, nu). Results are unchanged up to floating-point rounding.
* `spGLMexact()`, `spGLMstack()`, `recoverGLMscale()`: distances are computed from `coords` in C++, and `spGLMstack()` builds each correlation matrix once per (phi, nu).
* `spGLMexact()`, `spGLMstack()`, `stvcGLMexact()`, `stvcGLMstack()`: faster GLM backend. The pre-processing of the conjugate sampler needs one n x n Cholesky factorization instead of several, and is about 3x faster for `stvcGLMexact()`. Exact leave-one-out and K-fold cross-validation update the full-data pre-processing for each held-out site or fold (O(n^2) per site) instead of recomputing it (O(n^3)), with a direct recomputation as a fallback. Posterior and Monte Carlo draws are projected in blocks of 64 with level-3 BLAS, drawing the random variates in the original order: results are unchanged up to floating-point rounding, and speed-ups are largest with an optimized BLAS (e.g. about 3.5x for `spGLMexact()` and 5.5x for `stvcGLMexact()` with OpenBLAS). Results are unchanged up to floating-point rounding.
* `spGLMstack()`, `stvcGLMstack()`: candidate models that differ only in `boundary` share all the pre-processing (full data, and each leave-one-out site or cross-validation fold): candidates with the same (phi, nu), or (phi_s, phi_t), are fitted in one call that loops over their boundary values inside each fold. With three boundary values and 10-fold CV this cuts the run time by 25-55% with OpenBLAS (e.g. n = 1000: 23 s to 16 s for `spGLMstack()`, 26 s to 19 s for `stvcGLMstack()`; n = 2000: 146 s to 66 s for `spGLMstack()`) and by up to 15% with the reference BLAS. The posterior samples of these candidates are now drawn before their leave-one-out draws, so results for a given seed differ from earlier versions at the Monte Carlo level; a candidate fitted on its own (or a group with a single boundary value) is unchanged.
* `spGLMstack()`, `stvcGLMstack()`: new advanced tag `CV.update` in `loopd.controls` (`'auto'`, `'update'` or `'direct'`) choosing how the pre-processing of each cross-validation fold is obtained: by deletion updates of the full-data factors (scalar loops; faster with the reference BLAS) or by recomputing it on the fold (level-3 BLAS; faster with an optimized BLAS, e.g. a further 25-50% for `spGLMstack()` with OpenBLAS at n = 1000-2000). The default `'auto'` uses `'direct'` when R reports an optimized BLAS (OpenBLAS, MKL, BLIS, Accelerate/vecLib, ATLAS). Results agree up to floating-point rounding. `spGLMexact()` and `stvcGLMexact()` use `'auto'`.
* `spGLMstack()`, `stvcGLMstack()`: leaving out `loopd.controls`, or the `CV.K` or `nMC` tag (`spGLMstack()`), stopped with an error; the documented defaults are now used.
* `spGLMexact()`, `spGLMstack()`, `stvcGLMexact()`, `stvcGLMstack()`, `recoverGLMscale()`: a failed Cholesky factorization (e.g. a numerically singular correlation matrix from nearly coincident locations or very strong correlation, a prior covariance `V.beta` or `iw.scale` that is not positive definite) now stops with an informative error. Previously a message was printed to the console and the computation continued with invalid factors. In leave-one-out and cross-validation, a failed deletion update still falls back to the direct computation, and only a failure of that stops.
* `stvcGLMexact()` with `loopd = TRUE`, and `stvcGLMstack()`, for `process.type = "multivariate"`: the prior draw of the latent process in the posterior sampler used the Cholesky factor of the inverse-Wishart draw without zeroing its upper triangle, so the returned posterior samples of `beta`, `z` and `xi` had the wrong distribution. Fixed; the samples now match `stvcGLMexact()` with `loopd = FALSE` for the same seed. Leave-one-out predictive densities and stacking weights were not affected.
* `stvcGLMexact()`, `stvcGLMstack()`: K-fold cross-validation now saves R's random number state at the end, so later random numbers in the session no longer repeat the ones used by the cross-validation draws.
* `spLMexact()`, `spLMstack()`: new `diagnostics` element (one row per model) with the smallest relative Cholesky pivot of the n x n factorizations and the correlations of the two farthest-apart and the two closest locations, taken from quantities the fit computes anyway (no extra factorization). With `verbose = TRUE`, a "Diagnostics" section flags a nearly singular covariance matrix (pivot < 1e-8), an effective range far beyond the extent of the data (farthest-pair correlation > 0.95) or nearly uncorrelated locations (closest-pair correlation < 0.05); `spLMstack()` reports only flagged candidates with stacking weight above 0.05.
* GLMs: the latent pseudo-data are now drawn on the log scale, as log G1 - log G2 with gamma variates G1, G2 (binomial, binary) and as log G with an underflow-safe draw for shape < 1 (Poisson). Previously, a beta draw that rounded to 1, or a gamma draw that underflowed to 0, produced infinite or NaN samples and leave-one-out predictive densities; this happened with positive probability for small `boundary` (e.g. about 1 in 2 million draws at `boundary = 0.4` for binary data). The sampled distribution is unchanged, but results for a given seed differ from earlier versions at the Monte Carlo level (binomial and binary data; Poisson data with zero counts).
* `spGLMexact()`, `stvcGLMexact()`: `boundary` is no longer raised to 0.4 for binary data. A message is given when `boundary < 0.1`.
* `spGLMstack()`, `stvcGLMstack()`: each candidate `boundary` must lie in (0, 1), and a message is given when any is below 0.1.
* `stvcGLMexact()`, `stvcGLMstack()`: fixed a heap buffer overflow in K-fold cross-validation when the number of varying coefficients exceeded the number of fixed-effect covariates.
* `stvcGLMexact()`, `stvcGLMstack()`: `process.type` is now checked to be one of 'independent', 'independent.shared' or 'multivariate'. The undocumented, incomplete "multivariate2" option has been removed.
* `posteriorPredict()` and `recoverGLMscale()` for GLMs: faster and leaner (the kriging means of 64 draws at a time with level-3 BLAS, one n x n matrix fewer, the joint conditional covariance built only for joint prediction). The C++ code now saves and restores R's random number state, so `.Random.seed` advances and the draws no longer overlap with later random numbers in the session; pointwise prediction variances are clamped at 0 against rounding.
* `recoverGLMscale()` for `stvcGLMexact()`/`stvcGLMstack()` fits with `process.type = "independent"`: the posterior of each process's scale used the shape (nu.z + n*r)/2 of the shared-scale model; it is now (nu.z + n)/2. Recovered scales from earlier versions were too concentrated and biased downwards.
* `posteriorPredict()` for `stvcGLMexact()`/`stvcGLMstack()` fits with `process.type = "multivariate"` and `joint = FALSE`: `mu.pred` and `y.pred` were computed from the unscaled noise instead of the predicted `z.pred`; fixed.
* `spLMexact()`, `spLMstack()`: the posterior scale of the variance is computed in residual form, avoiding loss of precision when the signal is large relative to the noise; Cholesky failures now stop with an informative error.
* `cholUpdateDel()`, `cholUpdateDelBlock()`: indices are now checked to be single integers between 1 and `n`; an index of 0 previously returned a zero matrix or crashed R.
* All model-fitting functions now stop with an informative error if coordinates are duplicated; for `stvcGLMexact()` and `stvcGLMstack()`, a duplicate must coincide in both space and time.
* PSIS leave-one-out predictive densities (`loopd.method = "PSIS"`) are reimplemented in C++ following the loo package, with O(S) memory; Pareto k diagnostics are returned as `loopd.pareto_k`.
* `spLMexact()`, `spLMstack()`: the inverse-gamma prior is now placed on the measurement error variance `sigmaSq`; posterior samples of the spatial variance are returned as `sigmaSq.z`. The default prior is now `priors = "flat"`, i.e., p(beta, sigmaSq) proportional to 1/sigmaSq.

# spStack 1.1.3

* Documentation Update: Pan, Zhang, Bradley and Banerjee (2025) accepted at Bayesian Analysis.
* Migrate to CVXR 1.8.1 API and resolve CRAN Results errors.
* Add a fallback option to the optimization routine.

# spStack 1.1.2

* Documentation Update: Zhang, Tang and Banerjee (2025) accepted at JASA.

# spStack 1.1.1

* `lmulm_XTilde_VC()`, `lmulv_XTilde_VC()`: Fixed address sanitizer issue with string comparison with pointer to string literal.

# spStack 1.1.0

* `stvcGLMexact()`, `stvcGLMstack()`: New functions for spatially-temporally varying coefficients GLM.

* `posteriorPredict()`: New functions for posterior predictive inference using predictive stacking.

* `recoverGLMscale()`: Utility for recovering posterior samples of scale parameters in spatial and spatial-temporal GLMs.

# spStack 1.0.1

* Fix a memory leak issue.

# spStack 1.0.0

* Initial CRAN submission.
