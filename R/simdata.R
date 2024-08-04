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
#' @param n Sample size
#' @param beta \eqn{p}{p}-dimensional vector of fixed effects
#' @param cor.fn Character string corresponding a correlation function, valid
#' inputs are 'exponential' and 'matern'.
#' @param spParams Spatial process parameters
#' @param spvar Fixed value of spatial variance
#' @param deltasq Noise-to-spatial variance ratio
#' @importFrom stats dist runif
#' @export
sim_spData <- function(n, beta, cor.fn, spParams, spvar, deltasq) {

    S <- data.frame(s1 = runif(n, 0, 1), s2 = runif(n, 0, 1))
    D <- as.matrix(dist(S))
    V <- spvar * spCor(D, cor.fn, spParams)
    z <- rmvn(1, rep(0, n), V)
    nugget <- deltasq * spvar

    if (length(beta) == 1) {

        y <- beta + z + rnorm(n, mean = 0, sd = sqrt(nugget))
        dat <- cbind(S, y, z)
        names(dat) <- c("s1", "s2", "y", "z_true")

    } else if (length(beta) > 1) {

        p <- length(beta)
        X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
        y <- X %*% beta + z + rnorm(n, mean = 0, sd = sqrt(nugget))
        dat <- cbind(S, X[, -1], y, z)
        names(dat) = c("s1", "s2", paste("x", 1:(p - 1), sep = ""), "y",
                       "z_true")

    }

    return(dat)

}
