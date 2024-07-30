#' Univariate spatial linear mixed model
#'
#' @param formula a symbolic description of the regression model to be fit. See example below.
#' @param data an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spLM_stack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
#' @param cov.model a quoted keyword that specifies the correlation function used to model the spatial dependence structure among the observations. Supported covariance model key words are: \code{'exponential'} and \code{'matern'}. See below for details.
#' @param n.samples number of posterior samples to be generated
#' @param beta.prior.mean prior mean of beta
#' @param beta.prior.precision prior precision matrix of beta
#' @param phi spatial decay parameter
#' @param nu spatial smoothness parameter (Matern covariogram)
#' @param alpha nugget-to-spatial variance ratio
#' @param sigma.sq.prior.shape IG prior shape
#' @param sigma.sq.prior.rate IG prior rate
#' @param sp.effects spatial effects
#' @param verbose logical.
#' @param ... currently no additional arguments.
#' @export
bayesGeostatExact <- function(formula, data = parent.frame(), n.samples, beta.prior.mean, beta.prior.precision, coords, cov.model = "exponential",
    phi, nu, alpha, sigma.sq.prior.shape, sigma.sq.prior.rate, sp.effects = TRUE, verbose = TRUE, ...) {

    #################################################### Check for unused args
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for (i in elip.args) {
        if (!i %in% formal.args)
            warning("'", i, "' is not an argument")
    }

    if (missing(n.samples)) {
        stop("error: n.samples must be specified")
    }
    if (missing(beta.prior.mean)) {
        stop("error: beta.prior.mean must be specified")
    }
    if (missing(beta.prior.precision)) {
        stop("error: beta.prior.precision must be specified")
    }
    if (missing(coords)) {
        stop("error: coords must be specified")
    }
    if (missing(phi)) {
        stop("error: phi must be specified")
    }
    if (missing(alpha)) {
        stop("error: alpha must be specified")
    }
    if (missing(sigma.sq.prior.shape)) {
        stop("error: sigma.sq.prior.shape must be specified")
    }
    if (missing(sigma.sq.prior.rate)) {
        stop("error: sigma.sq.prior.rate must be specified")
    }

    if (!cov.model %in% c("gaussian", "exponential", "matern", "spherical")) {
        stop("error: specified cov.model '", cov.model, "' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }

    if (cov.model == "matern")
        if (missing(nu)) {
            stop("error: nu must be specified")
        }

    if (missing(formula)) {
        stop("error: formula must be specified")
    }
    ## if(class(formula) == 'formula'){
    if (inherits(formula, "formula")) {

        holder <- parseFormula(formula, data)
        Y <- holder[[1]]
        X <- as.matrix(holder[[2]])
        x.names <- holder[[3]]

    } else {
        stop("error: formula is misspecified")
    }

    p <- ncol(X)
    n <- nrow(X)

    ## check for dim problems
    if (length(beta.prior.mean) != p)
        stop(paste("error: beta.prior.mean must be of length ", p, sep = ""))

    beta.prior.precision <- as.matrix(beta.prior.precision)
    if (nrow(beta.prior.precision) != ncol(beta.prior.precision))
        stop("error: beta.prior.precision must be square")

    if (nrow(beta.prior.precision) != p)
        stop(paste("error: beta.prior.precision must be of dimension ", p, sep = ""))

    ## show summary of stuff
    if (verbose) {
        cat("-------------------------------------------------\n")
        cat("\tGeneral model description\n")
        cat("-------------------------------------------------\n")
        cat(paste("Model fit with ", n, " observations.\n", sep = ""))
        cat(paste("Number of covariates ", p, " (including intercept if specified).\n", sep = ""))
        cat(paste("Using the ", cov.model, " spatial correlation model.\n\n", sep = ""))
        cat("-------------------------------------------------\n")
        cat("\t\tSampling\n")
        cat("-------------------------------------------------\n")
    }



    D <- as.matrix(dist(coords))

    if (cov.model == "exponential") {
        R <- exp(-phi * D)
    } else if (cov.model == "matern") {
        R <- (D * phi)^nu/(2^(nu - 1) * gamma(nu)) * besselK(x = D * phi, nu = nu)
        diag(R) <- 1
    } else if (cov.model == "gaussian") {
        R <- exp(-1 * ((phi * D)^2))
    } else if (cov.model == "spherical") {
        R <- D
        R[TRUE] <- 1
        R[D > 0 & D < 1/phi] <- 1 - 1.5 * phi * D[D > 0 & D <= 1/phi] + 0.5 * ((phi * D[D > 0 & D <= 1/phi])^3)
        R[D >= 1/phi] <- 0
    } else {
        stop("error: specified cov.model '", cov.model, "' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }

    V.y <- R + alpha * diag(nrow(R))
    V.y.sqrt <- chol(V.y)
    V.y.inv <- chol2inv(V.y.sqrt)
    ttildeXtildeX <- t(X) %*% V.y.inv %*% X
    mu <- beta.prior.mean

    posterior.precision <- beta.prior.precision + ttildeXtildeX
    posterior.variance <- chol2inv(chol(posterior.precision))
    posterior.mean <- posterior.variance %*% (beta.prior.precision %*% mu + t(X) %*% V.y.inv %*% Y)

    sigma.sq.posterior.shape <- sigma.sq.prior.shape + n/2
    sigma.sq.posterior.rate <- sigma.sq.prior.rate + 0.5 * (t(mu) %*% beta.prior.precision %*% mu + t(Y) %*% V.y.inv %*% Y - t(posterior.mean) %*%
        posterior.precision %*% posterior.mean)

    return(list(posterior.mean, posterior.variance))
    # posterior.samples <- as.matrix(normalIGammaSampler(n.samples, posterior.mean, posterior.variance, sigma.sq.posterior.shape,
    # sigma.sq.posterior.rate))

    # beta.posterior.samples <- posterior.samples[,1:p] sigma.sq.posterior.samples <- posterior.samples[,(p+1)]
    # tau.sq.posterior.samples <- alpha*sigma.sq.posterior.samples posterior.samples <- cbind(beta.posterior.samples,
    # sigma.sq.posterior.samples, tau.sq.posterior.samples) colnames(posterior.samples) <- c(x.names,'sigma.sq', 'tau.sq')

    # cat(paste('Sampled: ',n.samples,' of ',n.samples,', ',100,'%\n', sep=''))

    # if(sp.effects){ cat('-------------------------------------------------\n') cat('\tRecovering spatial effects\n')
    # cat('-------------------------------------------------\n') w <- matrix(0, n, n.samples) R.inv <- chol2inv(chol(R)) V.sp <-
    # chol2inv(chol(R.inv + (1/alpha)*diag(nrow(R)))) resid.posterior <- matrix(rep(Y, times=n.samples), nrow=length(Y),
    # ncol=n.samples) - X%*%t(beta.posterior.samples) sp.posterior.mean <- (1/alpha)*t(V.sp%*%resid.posterior) V.sp.root <-
    # t(chol(V.sp)) ## chol returns 'upper-triangular'; so t(); for (s in 1:n.samples) { ## Using rmvnorm is slow - it calculates
    # the matrix square-root each time ## sp.effects[s,] <- rmvnorm(1, sp..posterior.mean[s,], sigma.sq.posterior.samples[s]*V.sp)
    # ## Instead use the pre-computed V.sp.root z <- rnorm(nrow(V.sp), 0, 1) w[,s] <- sp.posterior.mean[s,] +
    # sqrt(sigma.sq.posterior.samples[s])*V.sp.root%*%z } ## w <- matrix(0, n, n.samples) ## ## R.eigen <- eigen(R) ## R.vals.inv
    # <- 1/R.eigen$values ## ## R.vecs <- R.eigen$vectors ## ## sigma.sq <- posterior.samples[,'sigma.sq'] ## tau.sq <-
    # posterior.samples[,'tau.sq'] ## beta <- as.matrix(posterior.samples[,1:p]) ## ## R.vects.t <- t(R.vecs) ## ## for(s in
    # 1:n.samples){ ## ## S.w <- R.vecs%*%diag(1/(1/sigma.sq[s]*R.vals.inv+1/tau.sq[s]))%*%R.vects.t ## ## S.mu <- S.w%*%(Y -
    # X%*%as.matrix(beta[s,]))/tau.sq[s] ## ## S.w.sq <- R.vecs%*%diag(sqrt(1/(1/sigma.sq[s]*R.vals.inv+1/tau.sq[s]))) ## ## w[,s]
    # <- S.w.sq%*%as.matrix(rnorm(n))+S.mu ## } } ##end sp.effects ## ##return ## out <- list() out$X <- X out$n <- n out$p <- p
    # out$Y <- Y out$n.samples <- n.samples out$beta.prior.mean <- beta.prior.mean out$beta.prior.precision <- beta.prior.precision
    # out$coords <- coords out$cov.model <- cov.model out$phi <- phi out$alpha <- alpha out$sigma.sq.prior.shape <-
    # sigma.sq.prior.shape out$sigma.sq.prior.rate <- sigma.sq.prior.rate out$alpha <- alpha out$verbose <-verbose if(cov.model ==
    # 'matern') out$nu <- nu out$p.samples <- mcmc(posterior.samples) # if(sp.effects) # out$sp.effects <- w class(out) <-
    # 'bayesGeostatExact' out
}
