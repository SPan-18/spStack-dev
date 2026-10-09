# spStack (development version)

* `spLMexact()`, `spLMstack()`, `posteriorPredict()`: faster, leaner Gaussian backend. Exact leave-one-out predictive densities are computed in closed form, the fit holds a single n x n matrix, distances are computed from `coords` in C++, Matern correlations use closed forms for nu = 0.5, 1.5, 2.5, and `spLMstack()` builds each correlation matrix once per (phi, nu). Results are unchanged up to floating-point rounding.
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
