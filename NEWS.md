# spStack (development version)

* New model families: `svcLMexact()`, `svcLMstack()`, `svcGLMexact()` and `svcGLMstack()` for spatially varying coefficients model, and `stvcLMexact()`, `stvcLMstack()` for Gaussian spatially-temporally varying coefficients model.
* Faster, leaner and more robust backend: closed-form exact LOO for the Gaussian models, deletion updates for GLM LOO and K-fold CV, PSIS in C++, informative errors, and a new `diagnostics` element.
* Breaking changes: in `spLMexact()`/`spLMstack()` the inverse-gamma prior is on the noise variance `sigmaSq`; stacked fits store `model.params` and `stacking.summary`.
* New lean synthetic datasets `simSpatial` and `simSpaceTime` replace the old datasets and `sim_spData()`; all examples and vignettes use them. `surfaceplot()`/`surfaceplot2()` use a colorblind-friendly palette.

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
