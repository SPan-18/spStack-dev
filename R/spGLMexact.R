#' Univariate Bayesian spatial generalized linear model
#'
#' @description Fits a Bayesian spatial generalized linear model with fixed
#' values of spatial process parameters and some auxiliary model parameters. The
#' output contains posterior samples of the fixed effects, spatial random
#' effects and, if required, finds leave-one-out predictive densities.
#' @details With this function, we fit a Bayesian hierarchical spatial
#' generalized linear model by sampling exactly from the joint posterior
#' distribution utilizing the generalized conjugate multivariate distribution
#' theory (Bradley and Clinch 2024). Suppose \eqn{\chi = (s_1, \ldots, s_n)}
#' denotes the \eqn{n} spatial locations the response \eqn{y} is observed. Let
#' \eqn{y(s)} be the outcome at location \eqn{s} endowed with a probability law
#' from the natural exponential family, which we denote by
#' \deqn{
#' y(s) \sim \mathrm{EF}(x(s)^\top \beta + z(s); b, \psi)
#' }
#' for some positive parameter \eqn{b > 0} and unit log partition function
#' \eqn{\psi}. We consider the following response models based on the input
#' supplied to the argument `family`.
#' \describe{
#' \item{`'poisson'`}{It considers point-referenced Poisson responses
#' \eqn{y(s) \sim \mathrm{Poisson}(e^{x(s)^\top \beta + z(s)})}. Here,
#' \eqn{b = 1} and \eqn{\psi(t) = e^t}.}
#' \item{`'binomial'`}{It considers point-referenced binomial counts
#' \eqn{y(s) \sim \mathrm{Binomial}(m(s), \pi(s))} where, \eqn{m(s)} denotes
#' the total number of trials and probability of success
#' \eqn{\pi(s) = \mathrm{ilogit}(x(s)^\top \beta + z(s))} at location \eqn{s}.
#' Here, \eqn{b = m(s)} and \eqn{\psi(t) = \log(1+e^t)}.}
#' \item{`'binary'`}{It considers point-referenced binary data (0 or, 1) i.e.,
#' \eqn{y(s) \sim \mathrm{Bernoulli}(\pi(s))}, where probability of success
#' \eqn{\pi(s) = \mathrm{ilogit}(x(s)^\top \beta + z(s))} at location \eqn{s}.
#' Here, \eqn{b = 1} and \eqn{\psi(t) = \log(1 + e^t)}.}
#' }
#' The hierarchical model is given as
#' \deqn{
#' \begin{aligned}
#' y(s_i) &\mid \beta, z, \xi \sim EF(x(s_i)^\top \beta + z(s_i) +
#' \xi_i - \mu_i; b_i, \psi_y), i = 1, \ldots, n\\
#' \xi &\mid \beta, z, \sigma^2_\xi, \alpha_\epsilon \sim
#' \mathrm{GCM_c}(\cdots),\\
#' \beta &\mid \sigma^2_\beta \sim N(0, \sigma^2_\beta V_\beta), \quad
#' \sigma^2_\beta \sim \mathrm{IG}(\nu_\beta/2, \nu_\beta/2)\\
#' z &\mid \sigma^2_z \sim N(0, \sigma^2_z R(\chi; \phi, \nu)), \quad
#' \sigma^2_z \sim \mathrm{IG}(\nu_z/2, \nu_z/2),
#' \end{aligned}
#' }
#' where \eqn{\mu = (\mu_1, \ldots, \mu_n)^\top} denotes the discrepancy
#' parameter. We fix the spatial process parameters \eqn{\phi} and \eqn{\nu} and
#' the hyperparameters \eqn{V_\beta}, \eqn{\nu_\beta}, \eqn{\nu_z} and
#' \eqn{\sigma^2_\xi}. The term \eqn{\xi} is known as the fine-scale variation
#' term which is given a conditional generalized conjugate multivariate
#' distribution as prior. For more details, see Pan *et al.* 2024. Default
#' values for \eqn{V_\beta}, \eqn{\nu_\beta}, \eqn{\nu_z}, \eqn{\sigma^2_\xi}
#' are diagonal with each diagonal element 100, 2.1, 2.1 and 0.1 respectively.
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#'  If not found in \code{data}, the variables are taken from
#'  \code{environment(formula)}, typically the environment from which
#'  \code{spGLMexact} is called.
#' @param family Specifies the distribution of the response as a member of the
#'  exponential family. Supported options are `'poisson'`, `'binomial'` and
#'  `'binary'`.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors (optional) a list with each tag corresponding to a
#'  hyperparameter name and containing hyperprior details. Valid tags include
#'  `V.beta`, `nu.beta`, `nu.z` and `sigmaSq.xi`. Values of `nu.beta` and `nu.z`
#'  must be at least 2.1. If not supplied, uses defaults.
#' @param spParams fixed values of spatial process parameters.
#' @param boundary Specifies the boundary adjustment parameter. Must be a real
#' number between 0 and 1. Default is 0.5.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd logical. If `loopd=TRUE`, returns leave-one-out predictive
#'  densities, using method as given by \code{loopd.method}. Default is
#'  \code{FALSE}.
#' @param loopd.method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#'  valid inputs are `'exact'`, `'CV'` and `'PSIS'`. The option `'exact'`
#'  corresponds to exact leave-one-out predictive densities which requires
#'  computation almost equivalent to fitting the model \eqn{n} times. The
#'  options `'CV'` and `'PSIS'` are faster and they implement \eqn{K}-fold
#'  cross validation and Pareto-smoothed importance sampling to find approximate
#'  leave-one-out predictive densities (Vehtari *et al.* 2017).
#' @param CV.K An integer between 10 and 20. Considered only if
#' `loopd.method='CV'`. Default is 10 (as recommended in Vehtari *et. al* 2017).
#' @param loopd.nMC Number of Monte Carlo samples to be used to evaluate
#' leave-one-out predictive densities when `loopd.method` is set to either
#' 'exact' or 'CV'.
#' @param verbose logical. If \code{verbose = TRUE}, prints model description.
#' @param ... currently no additional argument.
#' @return An object of class \code{spGLMexact}, which is a list with the
#'  following tags -
#' \describe{
#' \item{priors}{details of the priors used, containing the values of the
#' boundary adjustment parameter (`boundary`), the variance parameter of the
#' fine-scale variation term (`simasq.xi`) and others.}
#' \item{samples}{a list of length 3, containing posterior samples of fixed
#'  effects (\code{beta}), spatial effects (\code{z}) and the fine-scale
#' variation term (\code{xi}).}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#'  densities.}
#' \item{model.params}{Values of the fixed parameters that includes
#'  \code{phi} (spatial decay), \code{nu} (spatial smoothness).}
#' }
#' The return object might include additional data that can be used for
#' subsequent prediction and/or model fit evaluation.
#' @seealso [spLMexact()]
#' @author Soumyakanti Pan <span18@ucla.edu>
#' @references Bradley JR, Clinch M (2024). "Generating Independent Replicates
#' Directly from the Posterior Distribution for a Class of Spatial Hierarchical
#' Models." *Journal of Computational and Graphical Statistics*, **0**(0), 1-17.
#' \doi{10.1080/10618600.2024.2365728}.
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2024). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' \doi{10.48550/arXiv.2406.04655}.
#' @references Vehtari A, Gelman A, Gabry J (2017). "Practical Bayesian Model
#' Evaluation Using Leave-One-out Cross-Validation and WAIC."
#' *Statistics and Computing*, **27**(5), 1413-1432. ISSN 0960-3174.
#' \doi{10.1007/s11222-016-9696-4}.
#' @examples
#' # Example 1: Analyze spatial poisson count data
#' data(simPoisson)
#' dat <- simPoisson[1:10, ]
#' mod1 <- spGLMexact(y ~ x1, data = dat, family = "poisson",
#'                    coords = as.matrix(dat[, c("s1", "s2")]),
#'                    cor.fn = "matern",
#'                    spParams = list(phi = 4, nu = 0.4),
#'                    n.samples = 100, verbose = TRUE)
#'
#' # summarize posterior samples
#' post_beta <- mod1$samples$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' # Example 2: Analyze spatial binomial count data
#' data(simBinom)
#' dat <- simBinom[1:10, ]
#' mod2 <- spGLMexact(cbind(y, n_trials) ~ x1, data = dat, family = "binomial",
#'                    coords = as.matrix(dat[, c("s1", "s2")]),
#'                    cor.fn = "matern",
#'                    spParams = list(phi = 3, nu = 0.4),
#'                    n.samples = 100, verbose = TRUE)
#'
#' # summarize posterior samples
#' post_beta <- mod2$samples$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' # Example 3: Analyze spatial binary data
#' data(simBinary)
#' dat <- simBinary[1:10, ]
#' mod3 <- spGLMexact(y ~ x1, data = dat, family = "binary",
#'                    coords = as.matrix(dat[, c("s1", "s2")]),
#'                    cor.fn = "matern",
#'                    spParams = list(phi = 4, nu = 0.4),
#'                    n.samples = 100, verbose = TRUE)
#'
#' # summarize posterior samples
#' post_beta <- mod3$samples$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#' @export
spGLMexact <- function(formula, data = parent.frame(), family,
                       coords, cor.fn, priors,
                       spParams, boundary = 0.5, n.samples,
                       loopd = FALSE, loopd.method = "exact", CV.K = 10,
                       loopd.nMC = 500, verbose = TRUE, ...){

  ##### check for unused args #####
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if (!i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

  ##### family #####
  if(missing(family)){
    stop("Family not specified")
  }else{
    if(!is.character(family)){
      stop("Family must be a character string. Choose from c('poisson',
           'binary', 'binomial').")
    }
    family <- tolower(family)
    if(!family %in% c('poisson', 'binary', 'binomial')){
      stop("Invalid family. Choose from c('poisson', 'binary', 'binomial').")
    }
  }

  ##### formula #####
  if(missing(formula)){
    stop("Formula must be specified")
  }
  if(inherits(formula, "formula")){
    holder <- parseFormula(formula, data)
    if(family == "binomial"){
      if(!dim(holder[[1L]])[2] == 2){
        stop("Response must be of the form cbind(y, n_trials).")
      }
      y <- as.numeric(holder[[1L]][, 1])
      n.binom <- as.numeric(holder[[1L]][, 2])
    }else{
      y <-  as.numeric(holder[[1L]])
    }
    X <- as.matrix(holder[[2]])
    X.names <- holder[[3]]
  } else {
    stop("Formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)

  if(family == "binary"){n.binom <- rep(1.0, n)}
  if(family == "poisson"){n.binom <- rep(0.0, n)}

  if(family == "poisson"){
    if(any(y < 0)){
      stop("family = 'poisson' but y contains negative values")
    }
    if(any(!floor(y) == y)){
      warning("family = 'poisson' but response contains non-integer values")
    }
  }else if(family == "binary"){
    if(any(y != 0 & y != 1)){
      stop("family = 'binary', response can only be either 0 or 1")
    }
  }else if(family == "binomial"){
    if(any(!floor(n.binom) == n.binom)){
      warning("family = 'binomial' but n_trials contains non-integer values")
    }
    if(any(y > n.binom)){
      stop("Number of successes exceeds number of trials in the data.")
    }
  }

  ## storage mode
  storage.mode(y) <- "double"
  storage.mode(n.binom) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ##### coords #####
  if(!is.matrix(coords)){
    stop("coords must n-by-2 matrix of xy-coordinate locations")
  }
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("either the coords have more than two columns or, number of rows is
         different than data used in the model formula")
  }

  coords.D <- 0
  coords.D <- iDist(coords)

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("cor.fn must be specified")
  }
  if(!cor.fn %in% c("exponential", "matern")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('exponential', 'matern')")
  }

  ##### priors #####
  nu.beta <- 0
  nu.z <- 0
  sigmaSq.xi <- 0
  missing.flag <- 0

  if(missing(priors)){
    V.beta <- diag(rep(100.0, p))
    nu.beta <- 2.1
    nu.z <- 2.1
    sigmaSq.xi <- 0.1
  }else{
    names(priors) <- tolower(names(priors))
    if(!'v.beta' %in% names(priors)){
      missing.flag <- missing.flag + 1
      V.beta <- diag(rep(100.0, p))
    }else{
      V.beta <- priors[["v.beta"]]
      if(!is.numeric(V.beta) || length(V.beta) != p^2){
        stop(paste("priors[['V.beta']] must be a ", p, "x", p,
                   " covariance matrix.", sep = ""))
      }
    }
    if(!'nu.beta' %in% names(priors)){
      missing.flag <- missing.flag + 1
      nu.beta <- 2.1
    }else{
      nu.beta <- priors[['nu.beta']]
      if(!is.numeric(nu.beta) || length(nu.beta) != 1){
        stop("priors[['nu.beta']] must be a single numeric value.")
      }
      if(nu.beta < 2.1){
        message("Supplied nu.beta is less than 2.1. Setting it to defaults.")
        nu.beta <- 2.1
      }
    }
    if(!'nu.z' %in% names(priors)){
      missing.flag <- missing.flag + 1
      nu.z <- 2.1
    }else{
      nu.z <- priors[['nu.z']]
      if(!is.numeric(nu.z) || length(nu.z) != 1){
        stop("priors[['nu.z']] must be a single numeric value.")
      }
      if(nu.z < 2.1){
        message("Supplied nu.z is less than 2.1. Setting it to defaults.")
        nu.z <- 2.1
      }
    }
    if(!'sigmasq.xi' %in% names(priors)){
      missing.flag <- missing.flag + 1
      sigmaSq.xi <- 0.1
    }else{
      sigmaSq.xi <- priors[['sigmasq.xi']]
      if(!is.numeric(sigmaSq.xi) || length(sigmaSq.xi) != 1){
        stop("priors[['sigmasq.xi']] must be a single numeric value.")
      }
      if(sigmaSq.xi < 0){
        stop("priors[['sigmaSq.xi']] must be positive real number.")
      }
    }
    if(missing.flag > 0){
    message("Some priors were not supplied. Using defaults.")
    }
  }

  ## storage mode
  storage.mode(nu.beta) <- "double"
  storage.mode(nu.z) <- "double"
  storage.mode(sigmaSq.xi) <- "double"

  ##### spatial process parameters #####
  phi <- 0
  nu <- 0

  if(missing(spParams)){
    stop("spParams (spatial process parameters) must be supplied.")
  }

  names(spParams) <- tolower(names(spParams))

  if(!"phi" %in% names(spParams)){
    stop("phi (spatial decay parameter) must be supplied.")
  }
  phi <- spParams[["phi"]]

  if(!is.numeric(phi) || length(phi) != 1){
    stop("phi (spatial decay parameter) must be a numeric scalar.")
  }
  if(phi <= 0){
    stop("phi (spatial decay parameter) must be a positive real number.")
  }

  if(cor.fn == "matern"){
    if (!"nu" %in% names(spParams)) {
      stop("nu (spatial smoothness parameter) must be supplied.")
    }
    nu <- spParams[["nu"]]
    if (!is.numeric(nu) || length(nu) != 1) {
      stop("nu (spatial smoothness parameter) must be a numeric scalar.")
    }
    if (nu <= 0) {
      stop("nu (spatial smoothness parameter) must be a positive real number.")
    }
  }

  ## storage mode
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"

  ##### Boundary adjustment parameter #####
  epsilon <- 0

  if(!is.numeric(boundary) || length(boundary) != 1){
    stop("boundary must be a numeric scalar between 0 and 1.")
  }else{
    epsilon <- boundary
    if(epsilon <= 0 || epsilon >= 1){
      message("boundary must be in the interval (0, 1). Using default of 0.5.")
      epsilon <- 0.5
    }
  }
  if(family == "binary"){
    if(epsilon < 0.4){
      message("family = binomial'; boundary < 0.4. Setting boundary = 0.4.")
      epsilon <- 0.4
    }
  }

  ## storage.mode
  storage.mode(epsilon) <- "double"

  ##### Leave-one-out setup #####

  if(loopd){
    if(missing(loopd.method)){
      stop("loopd.method must be specified")
    }else{
      loopd.method <- tolower(loopd.method)
    }
    if(!loopd.method %in% c("exact", "cv")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose from c('exact', 'CV').")
    }
    if(loopd.method == "exact"){
      CV.K <- as.integer(0)
    }
    if(loopd.method == "cv"){
      if(n < 100){
        message("Sample size too low for CV. Finding exact LOO-PD.")
        loopd.method <- "exact"
        CV.K <- as.integer(0)
      }else{
        if(CV.K < 10){
          message("CV.K must be at least 10. Setting it to 10.")
          CV.K <- 10
        }else if(CV.K > 20){
          message("CV.K must be at most 20. Setting it to 20.")
          CV.K <- 20
        }
        if(floor(CV.K) != CV.K){
          message("CV.K must be integer. Setting it to nearest integer.")
        }
      }
    }
    if(loopd.nMC < 500){
      message("Number of Monte Carlo samples too low. Using defaults = 500.")
      loopd.nMC = 500
    }
  }else{
    loopd.method <- "none"
  }

  storage.mode(CV.K) <- "integer"
  storage.mode(loopd.nMC) <- "integer"

  ##### sampling setup #####

  if(missing(n.samples)){
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ##### main function call #####
  ptm <- proc.time()

  if(loopd){
    samps <- .Call("spGLMexactLOO", y, X, p, n, family, n.binom, coords.D,
                   cor.fn, V.beta, nu.beta, nu.z, sigmaSq.xi, phi, nu, epsilon,
                   n.samples, loopd, loopd.method, CV.K, loopd.nMC, verbose)
  }else{
    samps <- .Call("spGLMexact", y, X, p, n, family, n.binom, coords.D, cor.fn,
                   V.beta, nu.beta, nu.z, sigmaSq.xi, phi, nu, epsilon,
                   n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$family <- family
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$priors <- list(mu.beta = rep(0, p), V.beta = V.beta, nu.beta = nu.beta,
                     nu.z = nu.z, sigmasq.xi = sigmaSq.xi, boundary = epsilon)
  out$n.samples <- n.samples
  out$samples <- samps[c("beta", "z", "xi")]
  if(loopd){
    out$loopd <- samps[["loopd"]]
    out$loopd.method <- loopd.method
    out$loopd.nMC <- loopd.nMC
  }
  if(cor.fn == 'matern'){
    out$model.params <- list(phi = phi, nu = nu)
  }else{
    out$model.params <- list(phi = phi)
  }
  out$run.time <- run.time

  class(out) <- "spGLMexact"

  return(out)

}