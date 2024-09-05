#' @description This package delivers Bayesian inference for point-referenced
#' Gaussian, poisson, binomial, and binary data using stacking of predictive
#' densities. Our algorithm is highly parallelisable and hence, much faster than
#' Markov chain Monte Carlo algorithms, and, delivers competitive predictive
#' performance. See Zhang, Tang, and Banerjee (2024), and, Pan, Zhang, Bradley,
#' and Banerjee (2024) for details.
#' 
#' The algorithm involves obtaining inference for a number of
#' candidate models specified by a grid of values of spatial process parameters
#' and some auxiliary model parameters. Inference of these individual models
#' proceed by sampling exactly from its joint posterior distribution in case of
#' both Gaussian and non-Gaussian data by utilizing the conjugate Bayesian
#' linear modelling framework and the generalized conjugate multivariate
#' distribution theory respectively.
#'
#' @details \tabular{ll}{ Package: \tab spStack\cr Type: \tab Package\cr
#' Version: \tab 0.1.0\cr Date: \tab 2024-09-03\cr License: \tab GPL-3\cr }
#' Accepts a formula, for example, \code{y~x1+x2}, for most regression models
#' accompanied by candidate values of spatial process parameters, and returns
#' posterior samples of the regression coefficients and the latent spatial
#' random effects. Posterior inference or prediction of any quantity of interest
#' will proceed from these samples. Main functions are - \cr [spLMexact()]\cr
#' [spLMstack()]\cr [spGLMexact()]\cr \code{spGLMstack()}
#'
#' @name spStack-package
#' @references Zhang L, Tang W, Banerjee S (2024). “Bayesian Geostatistics Using
#' Predictive Stacking.” \doi{10.48550/arXiv.2304.12414}.
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2024). “Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking.”
#' \doi{10.48550/arXiv.2406.04655}.
#' @keywords models spatial geostatistics Bayesian regression package
"_PACKAGE"