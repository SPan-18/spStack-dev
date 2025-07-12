#' Bayesian spatially-temporally varying coefficients generalized linear model
#' using predictive stacking
#'
#' @description Fits Bayesian spatial-temporal generalized linear model with
#' spatially-temporally varying coefficients on a collection of candidate models
#' constructed based on some candidate values of some model parameters specified
#' by the user and subsequently combines inference by stacking predictive
#' densities. See Pan, Zhang, Bradley, and Banerjee (2024) for more details.
#'
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially-temporally varying
#' coefficients. See examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{stvcGLMstack} is called.
#' @param family Specifies the distribution of the response as a member of the
#' exponential family. Supported options are `'poisson'`, `'binomial'` and
#' `'binary'`.
#' @param sp_coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' spatial coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param time_coords an \eqn{n \times 1}{n x 1} matrix of the observation
#' temporal coordinates in \eqn{\mathcal{T} \subseteq [0, \infty)}.
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial-temporal dependence structure among the observations.
#' Supported covariance model key words are: \code{'gneiting-decay'} (Gneiting
#' and Guttorp 2010). See below for details.
#' @param process.type a quoted keyword specifying the model for the
#' spatial-temporal process. Supported keywords are `'independent'` which
#' indicates independent processes for each varying coefficients characterized
#' by different process parameters, `independent.shared` implies independent
#' processes for the varying coefficients that shares common process parameters,
#' and `multivariate` implies correlated processes for the varying coefficients
#' modeled by a multivariate Gaussian process with an inverse-Wishart prior on
#' the correlation matrix. The input for `sptParams` and `priors` must be given
#' accordingly.
#' @param priors (optional) a list with each tag corresponding to a
#' hyperparameter name and containing hyperprior details. Valid tags include
#' `V.beta`, `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
#' and `nu.z` must be at least 2.1. If not supplied, uses defaults.
#' @param candidate.models an object of class `candidateModels` containing a
#' list of candidate models for stacking. See [candidateModels()] for details.
#' @param n.samples number of samples to be drawn from the posterior
#' distribution.
#' @param loopd.controls a list with details on how leave-one-out predictive
#' densities (LOO-PD) are to be calculated. Valid tags include `method`, `CV.K`
#' and `nMC`. The tag `method` can be either `'exact'` or `'CV'`. If sample size
#' is more than 100, then the default is `'CV'` with `CV.K` equal to its default
#' value 10 (Gelman *et al.* 2024). The tag `nMC` decides how many Monte Carlo
#' samples will be used to evaluate the leave-one-out predictive densities,
#' which must be at least 500 (default).
#' @param parallel logical. If \code{parallel=FALSE}, the parallelization plan,
#' if set up by the user, is ignored. If \code{parallel=TRUE}, the function
#' inherits the parallelization plan that is set by the user via the function
#' [future::plan()] only. Depending on the parallel backend available, users
#' may choose their own plan. More details are available at
#' \url{https://cran.R-project.org/package=future}.
#' @param solver (optional) Specifies the name of the solver that will be used
#' to obtain optimal stacking weights for each candidate model. Default is
#' \code{'ECOS'}. Users can use other solvers supported by the
#' \link[CVXR]{CVXR-package} package.
#' @param verbose logical. If \code{TRUE}, prints model-specific optimal
#' stacking weights.
#' @param ... currently no additional argument.
#' @return An object of class \code{stvcGLMstack}, which is a list including the
#'  following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 3, containing posterior
#'  samples of fixed effects (\code{beta}), spatial effects (\code{z}), and
#'  fine scale variation `xi` for that model.}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`candidate.models`}{a list of length \code{n_model} rows with each
#' entry containing details of the model parameters.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#' candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`solver.status`}{solver status as returned by the optimization
#' routine.}
#' }
#' This object can be further used to recover posterior samples of the scale
#' parameters in the model, and subsequrently, to make predictions at new
#' locations or times using the function [posteriorPredict()].
#' @importFrom rstudioapi isAvailable
#' @importFrom parallel detectCores
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_lapply
#' @examples
#' \donttest{
#' set.seed(1234)
#' data("sim_stvcPoisson")
#' dat <- sim_stvcPoisson[1:100, ]
#'
#' # create list of candidate models (multivariate)
#' mod.list2 <- candidateModels(list(phi_s = list(2, 3),
#'                                   phi_t = list(1, 2),
#'                                   boundary = c(0.5, 0.75)), "cartesian")
#'
#' # fit a spatial-temporal varying coefficient model using predictive stacking
#' mod1 <- stvcGLMstack(y ~ x1 + (x1), data = dat, family = "poisson",
#'                      sp_coords = as.matrix(dat[, c("s1", "s2")]),
#'                      time_coords = as.matrix(dat[, "t_coords"]),
#'                      cor.fn = "gneiting-decay",
#'                      process.type = "multivariate",
#'                      candidate.models = mod.list2,
#'                      loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
#'                      n.samples = 500)
#' }
#' @export
stvcGLMstack <- function(formula, data = parent.frame(), family,
                         sp_coords, time_coords, cor.fn, process.type, priors,
                         candidate.models, n.samples, loopd.controls,
                         parallel = FALSE, solver = "ECOS", verbose = TRUE, ...){

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
    holder <- parseFormula2(formula, data)
    if(family == "binomial"){
      if(!dim(holder[[1L]])[2] == 2){
        stop("Response must be of the form cbind(y, n_trials).")
      }
      y <- as.numeric(holder[[1L]][, 1])
      n.binom <- as.numeric(holder[[1L]][, 2])
    }else{
      y <-  as.numeric(holder[[1L]])
    }
    X <- as.matrix(holder[[2L]])
    X.names <- holder[[3L]]
    X_tilde <- as.matrix(holder[[5L]])
    X_tilde.names <- holder[[6L]]
    if(is.null(X_tilde)){
      stop("Formula does not indicate varying coefficients terms")
    }
  }else{
    stop("Formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(X_tilde)

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
  storage.mode(X_tilde) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(r) <- "integer"

  ##### coords #####
  if(!is.matrix(sp_coords)){
    stop("sp_coords must n-by-2 matrix of xy-coordinate locations")
  }
  if(ncol(sp_coords) != 2 || nrow(sp_coords) != n){
    stop("either the sp_coords have more than two columns or, number of rows is
         different than data used in the model formula")
  }

  if(!is.matrix(time_coords)){
    stop("time_coords must n-by-1 matrix of temporal coordinates")
  }
  if(ncol(time_coords) != 1 || nrow(time_coords) != n){
    stop("either the time_coords have more than one column or, number of rows is
         different than data used in the model formula")
  }

  storage.mode(sp_coords) <- "double"
  storage.mode(time_coords) <- "double"

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("cor.fn must be specified")
  }
  if(!cor.fn %in% c("gneiting-decay")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('gneiting-decay')")
  }

  ##### priors #####
  nu.beta <- 0
  nu.z <- 0
  sigmaSq.xi <- 0
  IW.scale <- 0
  missing.flag <- 0

  if(missing(priors)){
    V.beta <- diag(rep(100.0, p))
    nu.beta <- 2.1
    nu.z <- 2.1
    sigmaSq.xi <- 0.1
    IW.scale <- diag(rep(1.0, r))
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
    if(process.type == 'multivariate'){
      if(!'iw.scale' %in% names(priors)){
        missing.flag <- missing.flag + 1
        IW.scale <- diag(rep(1.0, r))
      }else{
        IW.scale <- priors[['iw.scale']]
        if(!is.numeric(IW.scale) || length(IW.scale) != r^2){
          stop(paste("priors[['IW.scale']] must be a ", r, "x", r,
                     " covariance matrix.", sep = ""))
        }
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

  #### set-up candidate.models for stacking parameters ####

  if(missing(candidate.models)){
    stop("error: candidate.models must be supplied.")
  }else{
    if(!inherits(candidate.models, "candidateModels")){
        stop("error: candidate.models must be an object of class 'candidateModels'.")
    }
    list_candidate <- candidate.models
    if(cor.fn == "gneiting-decay"){
      if(process.type %in% c("independent.shared", "multivariate")){
        check_validity <- all(vapply(list_candidate, function(x){
          length(x) == 3 &&
          identical(sort(names(x)), c("boundary", "phi_s", "phi_t")) &&
          all(vapply(x, function(vec) is.numeric(vec) && length(vec) == 1, logical(1)))
        }, logical(1)))
      }else if(process.type == "independent"){
        check_validity <- all(vapply(list_candidate, function(x){
          length(x) == 3 &&
          identical(sort(names(x)), c("boundary", "phi_s", "phi_t")) &&
          is.numeric(x$phi_s) && length(x$phi_s) == 2 &&
          is.numeric(x$phi_t) && length(x$phi_t) == 2 &&
          is.numeric(x$boundary) && length(x$boundary) == 1
          }, logical(1)))
      }
      if(!check_validity){
        stop("error: each element of candidate.models must be a list of length 3
             with names 'phi_s', 'phi_t' and 'boundary' containing numeric
             values. If process.type = 'independent', then 'phi_s' and 'phi_t'
             must be vectors of length 2, otherwise, they must be scalars.")
      }
    }
  }

  #### Leave-one-out setup ####
  loopd <- TRUE

  if(missing(loopd.controls)){
    if(n > 99){
      loopd.controls <- list()
      loopd.controls[["method"]] <- "CV"
      loopd.controls[["CV.K"]] <- 10
      loopd.controls[["nMC"]] <- 500
    }else{
      loopd.controls <- list()
      loopd.controls[["method"]] <- "exact"
      loopd.controls[["CV.K"]] <- 0
      loopd.controls[["nMC"]] <- 500
    }
  }else{
    names(loopd.controls) <- tolower(names(loopd.controls))
    if(!"method" %in% names(loopd.controls)){
      stop("error: method missing from loopd.controls.")
    }
    loopd.method <- loopd.controls[["method"]]
    loopd.method <- tolower(loopd.method)
    if(!loopd.method %in% c("exact", "cv")){
      stop("method = '", loopd.method, "' is not a valid option; choose from c('exact', 'CV').")
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
        if(!"cv.k" %in% names(loopd.controls)){
          message("CV.K missing from loopd.controls. Using defaults.")
          CV.K <- 10
        }
        CV.K <- loopd.controls[["cv.k"]]
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
    if(!"nmc" %in% names(loopd.controls)){
      message("nMC missing from loopd.controls. Using defaults.")
      loopd.nMC <- 500
    }else{
      loopd.nMC <- loopd.controls[["nmc"]]
    }
    if(loopd.nMC < 500){
      message("Number of Monte Carlo samples too low. Using defaults.")
      loopd.nMC = 500
    }
  }

  storage.mode(CV.K) <- "integer"
  storage.mode(loopd.nMC) <- "integer"

  ##### sampling setup #####

  if (missing(n.samples)) {
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  verbose_child <- FALSE
  storage.mode(verbose_child) <- "integer"

  #### main function call ####
  ptm <- proc.time()

  if(parallel){

    # Get current plan invoked by future::plan() by the user
    current_plan <- future::plan()
    NWORKERS.machine <- parallel::detectCores()
    NWORKERS.future <- future::nbrOfWorkers()

    if(NWORKERS.future >= NWORKERS.machine){
        stop(paste("error: Number of workers requested exceeds/matches machine
                   limit. Choose a value less than or equal to",
                   NWORKERS.machine - 1, "to avoid overcommitment of resources."
                   ))
    }

    if(rstudioapi::isAvailable()){
      # Check if the current plan is multicore
      if(inherits(current_plan, "multicore")){
        stop("\tThe 'multicore' plan is considered unstable when called from
        RStudio. Either run the script from terminal or, switch to a
        suitable plan, for example, 'multisession', 'cluster'. See
        https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
        for details.")
      }
    }else{
      if(.Platform$OS.type == "windows"){
        if(inherits(current_plan, "multicore")){
          stop("\t'multicore' is not supported by Windows due to OS limitations.
            Instead, use 'multisession' or, 'cluster' plan.")
        }
      }
    }

    samps <- future_lapply(1:length(list_candidate), function(x){
                    .Call("stvcGLMexactLOO", y, X, X_tilde, n, p, r, family,
                          n.binom, sp_coords, time_coords, cor.fn, V.beta,
                          nu.beta, nu.z, sigmaSq.xi, IW.scale, process.type,
                          as.numeric(list_candidate[[x]][["phi_s"]]),
                          as.numeric(list_candidate[[x]][["phi_t"]]),
                          as.numeric(list_candidate[[x]][["boundary"]]),
                          n.samples, loopd, loopd.method, CV.K, loopd.nMC,
                          verbose_child)}, future.seed = TRUE)

  }else{

    # Get current plan invoked by future::plan() by the user
    current_plan <- future::plan()
    if(!inherits(current_plan, "sequential")){
      message("Parallelization plan other than 'sequential' setup but parallel
      is set to FALSE. Ignoring parallelization plan.")
    }

    samps <- lapply(1:length(list_candidate), function(x){
                    .Call("stvcGLMexactLOO", y, X, X_tilde, n, p, r, family,
                          n.binom, sp_coords, time_coords, cor.fn, V.beta,
                          nu.beta, nu.z, sigmaSq.xi, IW.scale, process.type,
                          as.numeric(list_candidate[[x]][["phi_s"]]),
                          as.numeric(list_candidate[[x]][["phi_t"]]),
                          as.numeric(list_candidate[[x]][["boundary"]]),
                          n.samples, loopd, loopd.method, CV.K, loopd.nMC,
                          verbose_child)
                        })

  }

  loopd_mat <- do.call("cbind", lapply(samps, function(x) x[["loopd"]]))

  out_CVXR <- get_stacking_weights(loopd_mat, solver = solver)

  run.time <- proc.time() - ptm

  w_hat <- out_CVXR$weights
  w_hat <- as.numeric(w_hat)
  solver_status <- out_CVXR$status
  w_hat <- sapply(w_hat, function(x) max(0, x))
  w_hat <- w_hat / sum(w_hat)

  stack_out <- as.matrix(do.call("rbind", lapply(list_candidate, unlist)))
  stack_out <- cbind(stack_out, round(w_hat, 3))
  # colnames(stack_out) = c("phi_s", "phi_t", "boundary", "weight")
  colnames(stack_out)[dim(stack_out)[2]] = "weight"
  rownames(stack_out) = paste("Model", 1:nrow(stack_out))

  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", 1:length(list_candidate), sep = "")

  samps <- lapply(samps, function(x) x[c("beta", "z", "xi")])
  names(samps) <- paste("Model", 1:length(list_candidate), sep = "")

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$X.stvc.names <- X_tilde.names
  out$family <- family
  out$sp_coords <- sp_coords
  out$time_coords <- time_coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$n.samples <- n.samples
  out$candidate.models <- list_candidate
  out$priors <- list(mu.beta = rep(0, p), V.beta = V.beta, nu.beta = nu.beta,
                     nu.z = nu.z, sigmaSq.xi = sigmaSq.xi, IW.scale = IW.scale)
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.controls
  out$n.models <- length(list_candidate)
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$solver.status <- solver_status

  class(out) <- "stvcGLMstack"

  return(out)

}