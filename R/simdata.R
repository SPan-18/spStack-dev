#' @importFrom stats rnorm
rmvn <- function(n, mu = 0, V = matrix(1)) {

    p <- length(mu)
    if (any(is.na(match(dim(V), p))))
        stop("error: dimension mismatch.")
    D <- chol(V)
    t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))

}

#' Simulate spatial data on unit square
#'
#' @description Generates synthetic spatial data of different types where the
#' spatial co-ordinates are sampled uniformly on an unit square. Different types
#' include point-referenced Gaussian, Poisson, binomial and binary data. The
#' design includes an intercept and fixed covariates sampled from a standard
#' normal distribution.
#' @param n sample size.
#' @param beta a \eqn{p}{p}-dimensional vector of fixed effects.
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#' @param spParams a numeric vector containing spatial process parameters -
#' e.g., spatial decay and smoothness.
#' @param spvar value of spatial variance parameter.
#' @param deltasq value of noise-to-spatial variance ratio.
#' @param family a character specifying the distribution of the response as a
#' member of the exponential family. Valid inputs are `'gaussian'`, `'poisson'`,
#' `'binary'`, and `'binomial'`.
#' @param n_binom necessary only when `family = 'binomial'`. Must be a
#'  vector of length `n` that will specify the number of trials for each
#'  observation. If it is of length 1, then that value is considered to be the
#'  common value for the number of trials for all `n` observations.
#' @return a `data.frame` object containing the columns -
#' \describe{
#' \item{`s1, s2`}{2D-coordinates in unit square}
#' \item{`x1, x2, ...`}{covariates, not including intercept}
#' \item{`y`}{response}
#' \item{`n_trials`}{present only when binomial data is generated}
#' \item{`z_true`}{true spatial effects with which the data is generated}
#' }
#' @importFrom stats dist runif rpois rbinom
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @examples
#' set.seed(1729)
#' n <- 10
#' beta <- c(2, 5)
#' phi0 <- 2
#' nu0 <- 0.5
#' spParams <- c(phi0, nu0)
#' spvar <- 0.4
#' deltasq <- 1
#' sim1 <- sim_spData(n = n, beta = beta, cor.fn = "matern",
#'                    spParams = spParams, spvar = spvar, deltasq = deltasq,
#'                    family = "gaussian")
#' @export
sim_spData <- function(n, beta, cor.fn, spParams, spvar, deltasq, family,
                       n_binom){

  S <- data.frame(s1 = runif(n, 0, 1), s2 = runif(n, 0, 1))
  D <- as.matrix(dist(S))
  V <- spvar * spCor(D, cor.fn, spParams)
  z <- rmvn(1, rep(0, n), V)

  family <- tolower(family)

  if(family == "gaussian"){
    if(missing(deltasq)){
      stop("deltasq (noise-to-spatial variance ratio) not supplied.")
    }
    nugget <- deltasq * spvar
    if(length(beta) == 1){
      y <- beta + z + rnorm(n, mean = 0, sd = sqrt(nugget))
      dat <- cbind(S, y, z)
      names(dat) <- c("s1", "s2", "y", "z_true")
    }else{
      p <- length(beta)
      X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
      y <- X %*% beta + z + rnorm(n, mean = 0, sd = sqrt(nugget))
      dat <- cbind(S, X[, -1], y, z)
      names(dat) = c("s1", "s2", paste("x", 1:(p - 1), sep = ""), "y", "z_true")
    }
  }else if(family == "poisson"){
    if(length(beta) == 1){
      mu <- beta + z
      y <- sapply(1:n, function(x) rpois(n = 1, lambda = exp(mu[x])))
      dat <- cbind(S, y, z)
      names(dat) <- c("s1", "s2", "y", "z_true")
    }else{
      p <- length(beta)
      X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
      mu <- X %*% beta + z
      y <- sapply(1:n, function(x) rpois(n = 1, lambda = exp(mu[x])))
      dat <- cbind(S, X[, -1], y, z)
      names(dat) = c("s1", "s2", paste("x", 1:(p - 1), sep = ""), "y", "z_true")
    }
  }else if(family == "binomial"){
    if(missing(n_binom)){
      stop("error: n_binom must be specified.")
    }
    if(length(n_binom) == 1){
      binom_size <- rep(n_binom, n)
    }else if(length(n_binom) == n){
      binom_size <- n_binom
    }else{
      stop("error: n_binom must be a numeric vector of length 1 or ", n, ".\n")
    }
    if(length(beta) == 1){
      mu <- beta + z
      y <- sapply(1:n, function(x) rbinom(n = 1, size = binom_size[x],
                                          prob = ilogit(mu[x])))
      dat <- cbind(S, y, binom_size, z)
      names(dat) <- c("s1", "s2", "y", "n_trials", "z_true")
    }else{
      p <- length(beta)
      X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
      mu <- X %*% beta + z
      y <- sapply(1:n, function(x) rbinom(n = 1, size = binom_size[x],
                                          prob = ilogit(mu[x])))
      dat <- cbind(S, X[, -1], y, binom_size, z)
      names(dat) = c("s1", "s2", paste("x", 1:(p - 1), sep = ""), "y",
                     "n_trials", "z_true")
    }
  }else if(family == "binary"){
    if(length(beta) == 1){
      mu <- beta + z
      y <- sapply(1:n, function(x) rbinom(n = 1, size = 1,
                                          prob = ilogit(mu[x])))
      dat <- cbind(S, y, z)
      names(dat) <- c("s1", "s2", "y", "z_true")
    }else{
      p <- length(beta)
      X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
      mu <- X %*% beta + z
      y <- sapply(1:n, function(x) rbinom(n = 1, size = 1,
                                          prob = ilogit(mu[x])))
      dat <- cbind(S, X[, -1], y, z)
      names(dat) = c("s1", "s2", paste("x", 1:(p - 1), sep = ""), "y", "z_true")
    }
  }

  return(dat)

}
