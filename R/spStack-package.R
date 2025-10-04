#' @description
#' This package delivers functions to fit Bayesian hierarchical spatial process
#' models for point-referenced Gaussian, Poisson, binomial, and binary data
#' using stacking of predictive densities. It involves sampling from
#' analytically available posterior distributions conditional upon some
#' candidate values of the spatial process parameters for both Gaussian response
#' model as well as non-Gaussian responses, and, subsequently assimilate
#' inference from these individual posterior distributions using Bayesian
#' predictive stacking. Our algorithm is highly parallelizable and hence, much
#' faster than traditional Markov chain Monte Carlo algorithms while
#' delivering competitive predictive performance.
#'
#' In context of inference for spatial point-referenced data,
#' Bayesian hierarchical models involve latent spatial processes characterized
#' by spatial process parameters, which besides lacking substantive relevance in
#' scientific contexts, are also weakly identified and hence, impedes
#' convergence of MCMC algorithms. This motivates us to build methodology that
#' involves fast sampling from posterior distributions conditioned on a grid of
#' the weakly identified model parameters and combine the inference by stacking
#' of predictive densities (Yao *et. al* 2018). We exploit the Bayesian
#' conjugate linear modeling framework for the Gaussian case (Zhang, Tang and
#' Banerjee 2025) and the generalized conjugate multivariate distribution theory
#' (Pan, Zhang, Bradley and Banerjee 2025) to analytically derive the individual
#' posterior distributions.
#'
#' @details \tabular{ll}{ Package: \tab spStack\cr Type: \tab Package\cr
#' Version: \tab 1.1.0\cr License: \tab GPL-3\cr }
#' Accepts a formula, e.g., \code{y~x1+x2}, for most regression models
#' accompanied by candidate values of spatial process parameters, and returns
#' posterior samples of the regression coefficients and the latent spatial
#' random effects. Posterior inference or prediction of any quantity of interest
#' proceed from these samples. Main functions are - \cr [spLMexact()]\cr
#' [spGLMexact()]\cr [spLMstack()]\cr [spGLMstack()]
#'
#' @name spStack-package
#' @useDynLib spStack, .registration = TRUE
#' @references Zhang L, Tang W, Banerjee S (2025). "Bayesian Geostatistics Using
#' Predictive Stacking." *Journal of the American Statistical Association*,
#' **In press**. \doi{10.1080/01621459.2025.2566449}.
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' \doi{10.48550/arXiv.2406.04655}.
#' @references Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion)." *Bayesian
#' Analysis*, **13**(3), 917-1007. \doi{10.1214/17-BA1091}.
"_PACKAGE"