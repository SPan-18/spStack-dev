#' Univariate spatial linear mixed model using Bayesian predictive stacking
#'
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#'  If not found in \code{data}, the variables are taken from
#'  \code{environment(formula)}, typically the environment from which
#'  \code{spLMstack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors a list with each tag corresponding to a parameter name and
#'  containing prior details.
#' @param params_list a list containing candidate values of spatial process
#'  parameters and noise-to-spatial variance ratio. See example.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd_method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#'  valid inputs are `'exact'`, `'PSIS'`. The option `'exact'` corresponds to
#'  exact leave-one-out predictive densities which requires computation almost
#'  equivalent to fitting the model \eqn{n} times. The option `'PSIS'` is
#'  faster, as it finds approximate leave-one-out predictive densities using
#'  Pareto-smoothed importance sampling (Gelman *et al.* 2024).
#' @param parallel logical. If \code{FALSE}, the parallelization plan, if set up
#'  by the user is ignored. If \code{TRUE}, inherits the parallelization plan
#'  set by the user using the \code{future} package.
#' @param solver Specifies the name of the solver that will be used to obtain
#'  optimal stacking weights for each candidate model.
#' @param verbose logical. If \code{verbose = TRUE}, prints model-specific
#'  optimal stacking weights.
#' @param ... currently no additional argument.
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). “Pareto
#'  Smoothed Importance Sampling.” *Journal of Machine Learning Research*,
#'  **25**(72), 1–58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @importFrom rstudioapi isAvailable
#' @importFrom parallel detectCores
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_lapply
#' @export
spLMstack <- function(formula, data = parent.frame(), coords, cor.fn,
                      priors, params_list, n.samples, loopd_method,
                      parallel = FALSE, solver, verbose = TRUE, ...){

  ##### check for unused args #####
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for (i in elip.args) {
    if (!i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

  ##### formula #####
  if (missing(formula)) {
    stop("error: formula must be specified!")
  }

  if (inherits(formula, "formula")) {
    holder <- parseFormula(formula, data)
    y <- holder[[1L]]
    X <- as.matrix(holder[[2L]])
    X.names <- holder[[3L]]
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
  if (!is.matrix(coords)) {
    stop("error: coords must n-by-2 matrix of xy-coordinate locations")
  }
  if (ncol(coords) != 2 || nrow(coords) != n) {
    stop("error: either the coords have more than two columns or,
    number of rows is different than data used in the model formula")
  }

  coords.D <- 0
  coords.D <- iDist(coords)

  ##### correlation function #####
  if (missing(cor.fn)) {
    stop("error: cor.fn must be specified")
  }
  if (!cor.fn %in% c("exponential", "matern")) {
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
    c('exponential', 'matern').")
  }

  ##### priors #####
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.IG <- 0

  if(missing(priors)){

    warning("prior list not supplied, using defaults.")

  }else{

    names(priors) <- tolower(names(priors))

    ## Setup prior for beta
    if("beta.norm" %in% names(priors)){
      beta.Norm <- priors[["beta.norm"]]
      if(!is.list(beta.Norm) || length(beta.Norm) != 2){
        stop("error: beta.Norm must be a list of length 2")
      }
      if(length(beta.Norm[[1]]) != p){
        stop(paste("error: beta.Norm[[1]] must be a vector of length, ", p, ".",
                   sep = ""))
      }
      if(length(beta.Norm[[2]]) != p^2){
        stop(paste("error: beta.Norm[[2]] must be a ", p, "x", p,
                   " correlation matrix.", sep = ""))
      }
      beta.prior <- "normal"
    }

    ## Setup prior for sigma.sq
    if (!"sigma.sq.ig" %in% names(priors)) {
      stop("error: sigma.sq.IG must be specified")
    }
    sigma.sq.IG <- priors[["sigma.sq.ig"]]

    if (!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2) {
      stop("error: sigma.sq.IG must be a vector of length 2")
    }
    if (any(sigma.sq.IG <= 0)) {
      stop("error: sigma.sq.IG must be a positive vector of length 2")
    }

  }

  ## storage mode
  storage.mode(sigma.sq.IG) <- "double"

  #### set-up params_list for stacking parameters ####

  if(missing(params_list)){

    warning("warning: params_list must be supplied. Using defaults.")

    params_list <- vector(mode = "list", length = 3)
    names(params_list) <- c("phi", "nu", "noise_sp_ratio")
    ####### 1L REQUIRES AUTOMATION #######
    params_list[[1L]] <- c(3, 5, 10)
    params_list[[2L]] <- c(0.5, 1.0, 1.5)
    params_list[[3L]] <- c(0.25, 1.0, 2.0)

  }else{

    names(params_list) <- tolower(names(params_list))

    if(!"phi" %in% names(params_list)){
      stop("error: candidate values of phi must be specified in params_list.")
    }

    if(!"nu" %in% names(params_list)){
      warning("warning: candidate values of nu not specified. Using defaults
              c(0.5, 1, 1.5).")
      params_list[["nu"]] <- c(0.5, 1.0, 1.5)
    }

    if(!"noise_sp_ratio" %in% names(params_list)){
      warning("warning: candidate values of noise_sp_ratio not specified. Using
              defaults c(0.25, 1, 2).")
      params_list[["noise_sp_ratio"]] <- c(0.25, 1.0, 2.0)
    }

    params_list <- params_list[c("phi", "nu", "noise_sp_ratio")]

  }

  # setup parameters for candidate models based on cartesian product of
  # candidate values of each parameter
  list_candidate <- candidate_models(params_list)

  #### Leave-one-out setup ####
  loopd <- TRUE
  CV_K <- as.integer(0)
  storage.mode(CV_K) <- "integer"

  if(missing(loopd_method)){
    warning("loopd_method not specified. Using 'exact'.")
  }

  loopd_method <- tolower(loopd_method)

  if(!loopd_method %in% c("exact", "psis")){
    stop("error: Invalid loopd_method. Valid options are 'exact' and 'PSIS'.")
  }

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
                    .Call("spLMexactLOO", y, X, p, n, coords.D,
                          beta.prior, beta.Norm, sigma.sq.IG,
                          as.numeric(list_candidate[[x]]["phi"]),
                          as.numeric(list_candidate[[x]]["nu"]),
                          as.numeric(list_candidate[[x]]["noise_sp_ratio"]),
                          cor.fn, n.samples, loopd, loopd_method,
                          verbose_child)
                          }, future.seed = TRUE)

  }else{

    # Get current plan invoked by future::plan() by the user
    current_plan <- future::plan()
    if(!inherits(current_plan, "sequential")){
      warning("Parallelization plan other than 'sequential' setup but parallel
      is set to FALSE. Ignoring parallelization plan.")
    }

    samps <- lapply(1:length(list_candidate), function(x){
                    .Call("spLMexactLOO", y, X, p, n, coords.D,
                          beta.prior, beta.Norm, sigma.sq.IG,
                          as.numeric(list_candidate[[x]]["phi"]),
                          as.numeric(list_candidate[[x]]["nu"]),
                          as.numeric(list_candidate[[x]]["noise_sp_ratio"]),
                          cor.fn, n.samples, loopd, loopd_method,
                          verbose_child)
                        })

  }

  loopd_mat <- do.call("cbind", lapply(samps, function(x) x[["loopd"]]))
  loopd_mat[loopd_mat < -10] <- -10
#   return(loopd_mat)

  out_CVXR <- get_stacking_weights(loopd_mat, solver = solver)
#   out_CVXR <- stacking_weights(loopd_mat, solver = solver)
#   w_hat <- loo::stacking_weights(loopd_mat)
  run.time <- proc.time() - ptm

  w_hat <- out_CVXR$weights
  solver_status <- out_CVXR$status
  w_hat <- sapply(w_hat, function(x) max(0, x))
  w_hat <- w_hat / sum(w_hat)
#   w_hat <- as.numeric(w_hat)
#   solver_status <- "BFGS"

  if(verbose){
    stack_out <- as.matrix(do.call("rbind", lapply(list_candidate, unlist)))
    stack_out <- cbind(stack_out, round(w_hat, 3))
    colnames(stack_out) = c("phi", "nu", "noise_sp_ratio", "weight")
    rownames(stack_out) = paste("Model", 1:nrow(stack_out))
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$beta.prior.norm <- beta.Norm
  out$run.time <- run.time
  out$samples <- samps
  out$stacking_weights <- w_hat
  out$solver_status <- solver_status

  class(out) <- "spLMstack"

  return(out)

}