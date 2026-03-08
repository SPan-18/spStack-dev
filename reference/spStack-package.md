# spStack: Bayesian Geostatistics Using Predictive Stacking

This package delivers functions to fit Bayesian hierarchical spatial
process models for point-referenced Gaussian, Poisson, binomial, and
binary data using stacking of predictive densities. It involves sampling
from analytically available posterior distributions conditional upon
some candidate values of the spatial process parameters for both
Gaussian response model as well as non-Gaussian responses, and,
subsequently assimilate inference from these individual posterior
distributions using Bayesian predictive stacking. Our algorithm is
highly parallelizable and hence, much faster than traditional Markov
chain Monte Carlo algorithms while delivering competitive predictive
performance.

In context of inference for spatial point-referenced data, Bayesian
hierarchical models involve latent spatial processes characterized by
spatial process parameters, which besides lacking substantive relevance
in scientific contexts, are also weakly identified and hence, impedes
convergence of MCMC algorithms. This motivates us to build methodology
that involves fast sampling from posterior distributions conditioned on
a grid of the weakly identified model parameters and combine the
inference by stacking of predictive densities (Yao *et. al* 2018). We
exploit the Bayesian conjugate linear modeling framework for the
Gaussian case (Zhang, Tang and Banerjee 2025) and the generalized
conjugate multivariate distribution theory (Pan, Zhang, Bradley and
Banerjee 2025) to analytically derive the individual posterior
distributions.

## Details

|          |         |
|----------|---------|
| Package: | spStack |
| Type:    | Package |
| Version: | 1.1.0   |
| License: | GPL-3   |

Accepts a formula, e.g., `y~x1+x2`, for most regression models
accompanied by candidate values of spatial process parameters, and
returns posterior samples of the regression coefficients and the latent
spatial random effects. Posterior inference or prediction of any
quantity of interest proceed from these samples. Main functions are -  
[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)  
[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)  
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md)  
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)

## References

Zhang L, Tang W, Banerjee S (2025). "Bayesian Geostatistics Using
Predictive Stacking." *Journal of the American Statistical Association*,
**In press**.
[doi:10.1080/01621459.2025.2566449](https://doi.org/10.1080/01621459.2025.2566449)
.

Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian Inference for
Spatial-temporal Non-Gaussian Data Using Predictive Stacking." *Bayesian
Analysis*, **In Press**.
[doi:10.1214/25-BA1582](https://doi.org/10.1214/25-BA1582) .

Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to Average
Bayesian Predictive Distributions (with Discussion)." *Bayesian
Analysis*, **13**(3), 917-1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091) .

## See also

Useful links:

- <https://span-18.github.io/spStack-dev/>

- Report bugs at <https://github.com/SPan-18/spStack-dev/issues>

## Author

**Maintainer**: Soumyakanti Pan <span18@ucla.edu>
([ORCID](https://orcid.org/0009-0005-9889-7112))

Authors:

- Sudipto Banerjee <sudipto@ucla.edu>
