#' Univariate spatial generalized linear model
#'
#' @description Fits a Bayesian spatial generalized linear model with spatial
#'  process parameters and some auxiliary model parameters fixed to a value
#'  supplied by the user. The output contains posterior samples of the fixed
#'  effects, variance parameter, spatial random effects and, if required,
#'  leave-one-out predictive densities.
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#'  If not found in \code{data}, the variables are taken from
#'  \code{environment(formula)}, typically the environment from which
#'  \code{spGLMexact} is called.
#' @param family Specifies the distribution of the response as a member of the
#'  exponential family.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors (optional )a list with each tag corresponding to a
#'  hyperparameter name and containing hyperprior details. Tags include
#'  `V.beta`, `nu.beta`, `nu.z` and `sigmaSq.xi`. Default value for `V.beta` is
#'  diagonal with each entry 100. Defaults for `nu.beta`, `nu.z` and
#'  `sigmaSq.xi` are 2.1, 2.1 and 1 respectively. `nu.beta` and `nu.z` must be
#'  at least 2.1.
#' @param spParams fixed value of spatial process parameters.
#' @param boundary a real number greater than 0 and less than 1. Default is 0.5.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd logical. If `loopd=TRUE`, returns leave-one-out predictive
#'  densities, using method as given by \code{loopd.method}. Deafult is
#'  \code{FALSE}.
#' @param loopd.method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#'  valid inputs are `'exact'` and `'PSIS'`. The option `'exact'` corresponds to
#'  exact leave-one-out predictive densities which requires computation almost
#'  equivalent to fitting the model \eqn{n} times. The option `'PSIS'` is
#'  faster and finds approximate leave-one-out predictive densities using
#'  Pareto-smoothed importance sampling (Gelman *et al.* 2024).
#' @param verbose logical. If \code{verbose = TRUE}, prints model description.
#' @param ... currently no additional argument.
#' @seealso [spLMexact()]
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). “Pareto
#'  Smoothed Importance Sampling.” *Journal of Machine Learning Research*,
#'  **25**(72), 1–58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @export
spGLMexact <- function(formula, data = parent.frame(), family,
                       coords, cor.fn, priors,
                       spParams, boundary = 0.5, n.samples,
                       loopd = FALSE, loopd.method = "exact",
                       verbose = TRUE, ...){

  ##### check for unused args #####
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if (!i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

  ##### formula #####
  if(missing(formula)){
    stop("error: formula must be specified!")
  }

  if(inherits(formula, "formula")){
    holder <- parseFormula(formula, data)
    y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    X.names <- holder[[3]]
  } else {
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)

  ## storage mode
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ##### coords #####
  if(!is.matrix(coords)){
    stop("error: coords must n-by-2 matrix of xy-coordinate locations")
  }
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or,
    number of rows is different than data used in the model formula")
  }

  coords.D <- 0
  coords.D <- iDist(coords)

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("error: cor.fn must be specified")
  }
  if(!cor.fn %in% c("exponential", "matern")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('exponential', 'matern').")
  }

  ##### priors #####
  nu.beta <- 0
  nu.z <- 0
  sigmaSq.xi <- 0

  if(missing(priors)){

    V.beta <- diag(rep(100.0, p))
    nu.beta <- 2.1
    nu.z <- 2.1
    sigmaSq.xi <- 1.0

  }else{

    names(priors) <- tolower(names(priors))

    if(!'v.beta' %in% names(priors)){
      warning("V.beta not supplied. Using defaults.")
      V.beta <- diag(rep(100.0, p))
    }else{
      V.beta <- priors[["v.beta"]]
      if(!is.numeric(V.beta) || length(V.beta) != p^2){
        stop(paste("error: priors[['V.beta']] must be a ", p, "x", p,
                   " covariance matrix.", sep = ""))
      }
    }

    if(!'nu.beta' %in% names(priors)){
      warning("nu.beta not supplied. Using defaults.")
      nu.beta <- 2.1
    }else{
      nu.beta <- priors[['nu.beta']]
      if(!is.numeric(nu.beta) || length(nu.beta) != 1){
        stop("error: priors[['nu.beta']] must be a single numeric value.")
      }
      if(nu.beta < 2.1){
        warning("Supplied nu.beta is less than 2.1. Setting it to defaults.")
        nu.beta <- 2.1
      }
    }

    if(!'nu.z' %in% names(priors)){
      warning("nu.z not supplied. Using defaults.")
      nu.z <- 2.1
    }else{
      nu.z <- priors[['nu.z']]
      if(!is.numeric(nu.z) || length(nu.z) != 1){
        stop("error: priors[['nu.z']] must be a single numeric value.")
      }
      if(nu.z < 2.1){
        warning("Supplied nu.z is less than 2.1. Setting it to defaults.")
        nu.z <- 2.1
      }
    }

    if(!'sigmasq.xi' %in% names(priors)){
      warning("sigmaSq.xi not supplied. Using defaults.")
      sigmaSq.xi <- 1.0
    }else{
      sigmaSq.xi <- priors[['sigmasq.xi']]
      if(!is.numeric(sigmaSq.xi) || length(sigmaSq.xi) != 1){
        stop("error: priors[['sigmasq.xi']] must be a single numeric value.")
      }
      if(sigmaSq.xi < 0){
        stop("priors[['sigmaSq.xi']] must be positive real number.")
      }
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
    stop("phi must be supplied.")
  }
  phi <- spParams[["phi"]]

  if(!is.numeric(phi) || length(phi) != 1){
    stop("phi must be a numeric scalar.")
  }
  if(phi <= 0){
    stop("phi (decay parameter) must be a positive real number.")
  }

  if(cor.fn == "matern"){
    if (!"nu" %in% names(spParams)) {
      stop("nu (smoothness parameter) must be supplied.")
    }
    nu <- spParams[["nu"]]
    if (!is.numeric(nu) || length(nu) != 1) {
      stop("nu must be a numeric scalar.")
    }
    if (nu <= 0) {
      stop("nu (smoothness parameter) must be a positive real number.")
    }
  }

  ## storage mode
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"

  ##### Boundary adjustment parameter #####
  epsilon <- 0

  if(!is.numeric(boundary) || length(boundary) != 1){
    stop("error: boundary must be a single numeric value between 0 and 1.")
  }else{
    epsilon <- boundary
    if(epsilon <= 0 || epsilon >= 1){
      warning("boundary must be in the interval (0, 1). Using default of 0.5.")
      epsilon <- 0.5
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
    if(!loopd.method %in% c("exact", "psis")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose
           from c('exact', 'PSIS').")
    }
  }else{
    loopd.method <- "none"
  }

  ##### sampling setup #####

  if (missing(n.samples)) {
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ##### main function call #####
  ptm <- proc.time()

  if(loopd){
    samps <- list(beta = 0, z = 0)
  }else{
    samps <- .Call("spGLMexact", y, X, p, n, family, coords.D, cor.fn, V.beta,
                   nu.beta, nu.z, sigmaSq.xi, phi, nu, epsilon,
                   n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$beta.prior.norm <- list(rep(0, p), V.beta)
  out$samples <- samps[c("beta", "z")]
  if(loopd){
    out$loopd <- samps[["loopd"]]
  }
  if(cor.fn == 'matern'){
    out$model.params <- c(phi, nu)
  }else{
    out$model.params <- c(phi)
  }
  out$run.time <- run.time

  class(out) <- "spGLMexact"

  return(out)

}