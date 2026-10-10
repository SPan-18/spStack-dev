#' Bayesian spatial generalized linear model using predictive stacking
#'
#' @description Fits Bayesian spatial generalized linear model on a collection
#' of candidate models constructed based on some candidate values of some model
#' parameters specified by the user and subsequently combines inference by
#' stacking predictive densities. See Pan, Zhang, Bradley, and Banerjee (2025)
#' for more details.
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{spLMstack} is called.
#' @param family Specifies the distribution of the response as a member of the
#' exponential family. Supported options are `'poisson'`, `'binomial'` and
#' `'binary'`.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors (optional) a list with each tag corresponding to a parameter
#' name and containing prior details. Valid tags include `V.beta`, `nu.beta`,
#' `nu.z` and `sigmaSq.xi`.
#' @param candidate.models an object of class `candidateModels` containing a
#' list of candidate models for stacking. See [candidateModels()] for details.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd.controls a list with details on how leave-one-out predictive
#' densities (LOO-PD) are to be calculated. Valid tags include `method`, `CV.K`,
#' `nMC` and `CV.update`. The tag `method` can be either `'exact'` or `'CV'`. If
#' sample size is more than 100, then the default is `'CV'` with `CV.K` equal to
#' its default value 10 (Gelman *et al.* 2024). The tag `nMC` decides how many
#' Monte Carlo samples will be used to evaluate the leave-one-out predictive
#' densities, which must be at least 500 (default). The tag `CV.update` is an
#' advanced option, used only if `method = 'CV'`, that decides how the
#' pre-processing of the model fit on each fold is obtained, and should be
#' changed with care as the faster choice depends on the BLAS library that R is
#' linked with (see `sessionInfo()`). `CV.update = 'update'` obtains it by
#' deletion updates of the full-data Cholesky factors, which is the faster
#' choice with the reference BLAS that R ships with. `CV.update = 'direct'`
#' recomputes it on each fold, which is faster only with an optimized BLAS such
#' as OpenBLAS, Intel MKL or Apple Accelerate (vecLib), and is slower
#' otherwise. The default `CV.update = 'auto'` uses `'direct'` if such an
#' optimized BLAS is detected from the library paths reported by R and
#' `'update'` otherwise; set it explicitly if the BLAS is not detected correctly
#' (for example, an optimized BLAS installed in place of `Rblas.dll` on
#' Windows). Both choices give the same results up to floating-point rounding,
#' so only the run time is affected.
#' @param parallel logical. If \code{parallel=FALSE}, the parallelization plan,
#'  if set up by the user, is ignored. If \code{parallel=TRUE}, the function
#'  inherits the parallelization plan that is set by the user via the function
#'  [future::plan()] only. Depending on the parallel backend available, users
#'  may choose their own plan. More details are available at
#'  \url{https://cran.R-project.org/package=future}.
#' @param solver (optional) Specifies the name of the solver that will be used
#'  to obtain optimal stacking weights for each candidate model. Default order
#'  is \code{c("CLARABEL", "ECOS", "SCS")}. Users can use other solvers
#'  supported by the \link[CVXR]{CVXR-package} package.
#' @param verbose logical. If \code{TRUE}, prints model-specific optimal
#' stacking weights.
#' @param ... currently no additional argument.
#' @return An object of class \code{spGLMstack}, which is a list including the
#'  following tags -
#' \describe{
#' \item{`family`}{the distribution of the responses as indicated in the
#' function call}
#' \item{`samples`}{a list of length equal to total number of candidate models
#' with each entry corresponding to a list of length 3, containing posterior
#' samples of fixed effects (\code{beta}), spatial effects (\code{z}) and
#' fine-scale variation term (\code{xi}) for that particular model.}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`loopd.method`}{a list containing details of the algorithm used for
#' calculation of leave-one-out predictive densities. For K-fold
#' cross-validation, its tag `cv.update` records the pre-processing method
#' (`'update'` or `'direct'`) that was used.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`model.params`}{a list with one element per candidate model, each a
#'  named list of its parameters: \code{phi}, \code{nu} (\code{NA} for the
#'  exponential correlation function) and \code{boundary} (boundary adjustment parameter).}
#' \item{`stacking.summary`}{a matrix with one row per candidate model,
#'  containing its parameters and its optimal stacking weight, for display.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#'  candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`diagnostics`}{a list of diagnostics. Element \code{numerical} is a
#' data frame with one row per candidate model and columns \code{min.pivot}
#' (the smallest relative Cholesky pivot of the \eqn{n \times n}{n x n}
#' correlation matrix; values below 1e-8 indicate a nearly singular
#' correlation matrix), \code{min.cor} and \code{max.cor} (the correlations
#' of the two farthest-apart and of the two closest locations; values of
#' \code{min.cor} above 0.95 suggest an effective range far exceeding the
#' extent of the data, values of \code{max.cor} below 0.05 nearly
#' uncorrelated locations), obtained from quantities the fit computes
#' anyway. Element \code{solver} describes the optimization for the stacking
#' weights: the solver used (\code{used}) and its status (\code{status}),
#' the installed and requested solvers, the search order, the attempts with
#' their status, and whether the fallback \code{loo::stacking_weights()} was
#' used. If \code{verbose = TRUE}, a "Diagnostics" section is printed if
#' there is an issue: numerical flags of candidate models with stacking
#' weight above 0.05 (extreme candidates with negligible weight are expected
#' in a stacking grid and are only counted), and solver problems (a
#' requested solver not installed, an inaccurate solution, or the
#' fallback).}
#' }
#' The return object might include additional data that is useful for subsequent
#' prediction, model fit evaluation and other utilities.
#' @details Instead of assigning a prior on the process parameters \eqn{\phi}
#' and \eqn{\nu}, the boundary adjustment parameter \eqn{\epsilon}, we consider
#' a set of candidate models based on some candidate values of these parameters
#' supplied by the user. Suppose the set of candidate models is
#' \eqn{\mathcal{M} = \{M_1, \ldots, M_G\}}. Then for each
#' \eqn{g = 1, \ldots, G}, we sample from the posterior distribution
#' \eqn{p(\sigma^2, \beta, z \mid y, M_g)} under the model \eqn{M_g} and find
#' leave-one-out predictive densities \eqn{p(y_i \mid y_{-i}, M_g)}. Then we
#' solve the optimization problem
#' \deqn{
#' \begin{aligned}
#' \max_{w_1, \ldots, w_G}& \, \frac{1}{n} \sum_{i = 1}^n \log \sum_{g = 1}^G
#' w_g p(y_i \mid y_{-i}, M_g) \\
#' \text{subject to} & \quad w_g \geq 0, \sum_{g = 1}^G w_g = 1
#' \end{aligned}
#' }
#' to find the optimal stacking weights \eqn{\hat{w}_1, \ldots, \hat{w}_G}.
#' @seealso [spGLMexact()], [spLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' *Bayesian Analysis*, **In Press**. \doi{10.1214/25-BA1582}.
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto
#'  Smoothed Importance Sampling." *Journal of Machine Learning Research*,
#'  **25**(72), 1-58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @importFrom rstudioapi isAvailable
#' @importFrom parallel detectCores
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_lapply
#' @examples
#' \donttest{
#' set.seed(1234)
#' data(simSpatial)
#' dat <- simSpatial[1:100, ]
#' cand.mod <- candidateModels(list(phi = c(3, 6, 10), nu = c(0.5, 1),
#'                                  boundary = c(0.5, 0.6)), "cartesian")
#'
#' mod1 <- spGLMstack(y_pois ~ x1 + x2, data = dat, family = "poisson",
#'                    coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
#'                    candidate.models = cand.mod,
#'                    n.samples = 1000,
#'                    loopd.controls = list(method = "CV", CV.K = 10, nMC = 1000),
#'                    parallel = TRUE, verbose = TRUE)
#'
#' post_samps <- stackedSampler(mod1)
#' post_beta <- post_samps$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' # compare the posterior medians of the spatial effects with the truth
#' cor(apply(post_samps$z, 1, median), dat$z_true)
#' }
#' @export
spGLMstack <- function(formula, data = parent.frame(), family,
                       coords, cor.fn, priors,
                       candidate.models, n.samples, loopd.controls,
                       parallel = FALSE, solver = NULL, verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

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

  check_no_missing(y = y, X = X, n.binom = n.binom, coords = coords)
  check_distinct_coords(coords)

  ## distances are computed in C++ from the coordinates
  storage.mode(coords) <- "double"

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
      check_cov_matrix(V.beta, p, "priors[[\'V.beta\']]")
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

  #### set-up candidate.models for stacking parameters ####

  if(missing(candidate.models)){
    stop("error: candidate.models must be supplied.")
  }else{
    if(!inherits(candidate.models, "candidateModels")){
      stop("error: candidate.models must be an object of class 'candidateModels'.")
    }
    if(cor.fn == "matern"){
      check_validity <- all(vapply(candidate.models, function(x){
        length(x) == 3 &&
        identical(sort(names(x)), c("boundary", "nu", "phi")) &&
        all(vapply(x, function(v) is.numeric(v) && length(v) == 1, logical(1)))
      }, logical(1)))
      if(!check_validity){
        stop("error: each element of candidate.models must be a named list with
             scalar numeric entries 'phi', 'nu' and 'boundary'.")
      }
    }else{
      check_validity <- all(vapply(candidate.models, function(x){
        length(x) == 2 &&
        identical(sort(names(x)), c("boundary", "phi")) &&
        all(vapply(x, function(v) is.numeric(v) && length(v) == 1, logical(1)))
      }, logical(1)))
      if(!check_validity){
        stop("error: each element of candidate.models must be a named list with
             scalar numeric entries 'phi' and 'boundary'.")
      }
      candidate.models <- lapply(candidate.models, function(x){
        x[["nu"]] <- NA_real_                                  # not used by the exponential
        x
      })
      class(candidate.models) <- "candidateModels"
    }
  }

  # boundary adjustment parameters of the candidate models
  cand_eps <- vapply(candidate.models, function(x) as.numeric(x[["boundary"]]), numeric(1))
  if(any(!is.finite(cand_eps) | cand_eps <= 0 | cand_eps >= 1)){
    stop("error: each 'boundary' in candidate.models must be in the interval (0, 1).")
  }
  if(any(cand_eps < 0.1)){
    message("candidate.models contains boundary < 0.1: the latent pseudo-data of observations at the edge of the support (y = 0, or y = trials) are very diffuse (standard deviation about 1/boundary on the linear predictor scale). Larger boundary values are recommended.")
  }

  list_candidate <- candidate.models

  #### Leave-one-out setup ####
  loopd <- TRUE

  # defaults if loopd.controls is not supplied; parsed below like a user-supplied list
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
  }

  if(!is.list(loopd.controls)){
    stop("error: loopd.controls must be a list.")
  }
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
      }else{
        CV.K <- loopd.controls[["cv.k"]]
      }
      if(CV.K < 10){
        message("CV.K must be at least 10. Setting it to 10.")
        CV.K <- 10
      }else if(CV.K > 20){
        message("CV.K must be at most 20. Setting it to 20.")
        CV.K <- 20
      }
      if(floor(CV.K) != CV.K){
        message("CV.K must be integer. Setting it to nearest integer.")
        CV.K <- round(CV.K)
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
    message("Number of Monte Carlo samples too low. Using defaults = 500.")
    loopd.nMC = 500
  }

  # pre-processing of the K-fold CV subsets: deletion updates or direct recomputation (advanced)
  if(!"cv.update" %in% names(loopd.controls)){
    CV.update <- "auto"
  }else{
    CV.update <- loopd.controls[["cv.update"]]
  }
  CV.update <- resolve_CV_update(CV.update)
  if(loopd.method == "cv"){
    loopd.controls[["cv.update"]] <- CV.update
  }
  CV.update <- as.integer(CV.update == "update")

  storage.mode(CV.K) <- "integer"
  storage.mode(loopd.nMC) <- "integer"

  ##### sampling setup #####

  if (missing(n.samples)) {
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  # candidate models sharing (phi, nu) are fitted in one call, which builds the
  # spatial correlation matrix and all the pre-processing (full data and
  # leave-one-out/CV subsets) once and shares it across their boundary values
  cand_phi <- vapply(list_candidate, function(x) as.numeric(x[["phi"]]), numeric(1))
  cand_nu <- vapply(list_candidate, function(x) as.numeric(x[["nu"]]), numeric(1))
  cand_boundary <- vapply(list_candidate, function(x) as.numeric(x[["boundary"]]),
                          numeric(1))
  cand_key <- sprintf("%a %a", cand_phi, cand_nu)          # exact (hexadecimal) keys
  cand_groups <- split(seq_along(list_candidate),
                       factor(cand_key, levels = unique(cand_key)))
  names(cand_groups) <- NULL

  # fits the candidate models in group g
  fit_group <- function(g){
    idx <- cand_groups[[g]]
    .Call(C_spGLMexactLOOgrid, y, X, p, n, family, n.binom,
          coords, cor.fn, V.beta, nu.beta, nu.z, sigmaSq.xi,
          cand_phi[idx[1]], cand_nu[idx[1]], cand_boundary[idx],
          n.samples, loopd, loopd.method, CV.K, loopd.nMC, CV.update)
  }

  # results of the groups, put back in the order of list_candidate
  ungroup <- function(samps_grouped){
    samps <- vector("list", length(list_candidate))
    for(g in seq_along(cand_groups)){
      samps[cand_groups[[g]]] <- samps_grouped[[g]]
    }
    samps
  }

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

    samps <- ungroup(future_lapply(seq_along(cand_groups), fit_group,
                                   future.seed = TRUE))

  }else{

    # Get current plan invoked by future::plan() by the user
    current_plan <- future::plan()
    if(!inherits(current_plan, "sequential")){
      message("Parallelization plan other than 'sequential' setup but parallel
      is set to FALSE. Ignoring parallelization plan.")
    }

    samps <- ungroup(lapply(seq_along(cand_groups), fit_group))

  }

  loopd_mat <- do.call("cbind", lapply(samps, function(x) x[["loopd"]]))

  # the solver details are kept in the 'diagnostics' element (and reported there
  # if there is an issue) instead of being printed
  out <- suppressMessages(get_stacking_weights(
    loopd_mat,
    solver = solver,
    verbose = FALSE
  ))

  w_hat <- out$weights
  if(identical(out$solver, "none")){
    message(loo_fallback_message())
  }
  solver_status <- out$status
  solver_used <- out$solver
  out_solver_details <- out$details

  run.time <- proc.time() - ptm

  # columns in a fixed order (exponential candidates carry nu = 0, appended last)
  stack_out <- as.matrix(do.call("rbind", lapply(list_candidate, function(x){
    unlist(x[c("phi", "nu", "boundary")])
  })))
  stack_out <- cbind(stack_out, round(w_hat, 3))
  colnames(stack_out) = c("phi", "nu", "boundary", "weight")
  rownames(stack_out) = paste("Model", 1:nrow(stack_out))

  if(cor.fn == 'exponential'){
    stack_out <- stack_out[, c("phi", "boundary", "weight")]
  }

  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", 1:length(list_candidate), sep = "")

  diagnostics <- list(numerical = collect_diagnostics(samps, paste("Model", seq_along(list_candidate))))
  diagnostics$solver <- c(list(used = solver_used, status = solver_status), out_solver_details)
  if(verbose){
    print_diagnostics(diagnostics, weights = w_hat)
  }

  samps <- lapply(samps, function(x) x[c("beta", "z", "xi")])
  names(samps) <- paste("Model", 1:length(list_candidate), sep = "")

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$family <- family
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$priors <- list(mu.beta = rep(0, p), V.beta = V.beta, nu.beta = nu.beta,
                     nu.z = nu.z, sigmasq.xi = sigmaSq.xi)
  out$n.samples <- n.samples
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.controls
  out$n.models <- length(list_candidate)
  # model parameters of each candidate, read by name (e.g. by posteriorPredict())
  out$model.params <- lapply(list_candidate, function(x){
    list(phi = as.numeric(x[["phi"]]), nu = as.numeric(x[["nu"]]),
         boundary = as.numeric(x[["boundary"]]))
  })
  names(out$model.params) <- paste("Model", seq_along(list_candidate))
  # table of the candidate models and their stacking weights (for display)
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$diagnostics <- diagnostics

  class(out) <- "spGLMstack"

  return(out)

}