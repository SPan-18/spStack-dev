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
#' @param params.list a list containing candidate values of spatial process
#' parameters for the `cor.fn` used, and, the boundary parameter.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd.controls a list with details on how leave-one-out predictive
#' densities (LOO-PD) are to be calculated. Valid tags include `method`, `CV.K`
#' and `nMC`. The tag `method` can be either `'exact'` or `'CV'`. If sample size
#' is more than 100, then the default is `'CV'` with `CV.K` equal to its default
#' value 10 (Gelman *et al.* 2024). The tag `nMC` decides how many Monte Carlo
#' samples will be used to evaluate the leave-one-out predictive densities,
#' which must be at least 500 (default).
#' @param parallel logical. If \code{parallel=FALSE}, the parallelization plan,
#'  if set up by the user, is ignored. If \code{parallel=TRUE}, the function
#'  inherits the parallelization plan that is set by the user via the function
#'  [future::plan()] only. Depending on the parallel backend available, users
#'  may choose their own plan. More details are available at
#'  \url{https://cran.R-project.org/package=future}.
#' @param solver (optional) Specifies the name of the solver that will be used
#' to obtain optimal stacking weights for each candidate model. Default is
#' \code{'ECOS'}. Users can use other solvers supported by the
#' \link[CVXR]{CVXR-package} package.
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
#' calculation of leave-one-out predictive densities.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`candidate.models`}{a matrix with \code{n_model} rows with each row
#'  containing details of the model parameters and its optimal weight.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#'  candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`solver.status`}{solver status as returned by the optimization
#' routine.}
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
#' \doi{10.48550/arXiv.2406.04655}.
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
#' data("simPoisson")
#' dat <- simPoisson[1:100,]
#' mod1 <- spGLMstack(y ~ x1, data = dat, family = "poisson",
#'                    coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
#'                   params.list = list(phi = c(3, 7, 10), nu = c(0.25, 0.5, 1.5),
#'                                      boundary = c(0.5, 0.6)),
#'                   n.samples = 1000,
#'                   loopd.controls = list(method = "CV", CV.K = 10, nMC = 1000),
#'                   parallel = TRUE, solver = "ECOS", verbose = TRUE)
#'
#' # print(mod1$solver.status)
#' # print(mod1$run.time)
#'
#' post_samps <- stackedSampler(mod1)
#' post_beta <- post_samps$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' post_z <- post_samps$z
#' post_z_summ <- t(apply(post_z, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
#'
#' z_combn <- data.frame(z = dat$z_true,
#'                       zL = post_z_summ[, 1],
#'                       zM = post_z_summ[, 2],
#'                       zU = post_z_summ[, 3])
#'
#' library(ggplot2)
#' plot_z <- ggplot(data = z_combn, aes(x = z)) +
#'  geom_errorbar(aes(ymin = zL, ymax = zU),
#'                width = 0.05, alpha = 0.15,
#'                color = "skyblue") +
#'  geom_point(aes(y = zM), size = 0.25,
#'             color = "darkblue", alpha = 0.5) +
#'  geom_abline(slope = 1, intercept = 0,
#'              color = "red", linetype = "solid") +
#'  xlab("True z") + ylab("Posterior of z") +
#'  theme_bw() +
#'  theme(panel.background = element_blank(),
#'        aspect.ratio = 1)
#' }
#' @export
spGLMstack <- function(formula, data = parent.frame(), family,
                       coords, cor.fn, priors,
                       params.list, n.samples, loopd.controls,
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

  #### set-up params.list for stacking parameters ####

  if(missing(params.list)){

    stop("error: params.list must be supplied.")

  }else{

    names(params.list) <- tolower(names(params.list))

    if(!"phi" %in% names(params.list)){
      stop("error: candidate values of phi must be specified in params_list.")
    }

    if(cor.fn == 'matern'){
      if(!"nu" %in% names(params.list)){
        message("Candidate values of nu not specified. Using defaults
                c(0.5, 1, 1.5).")
        params.list[["nu"]] <- c(0.5, 1.0, 1.5)
      }
    }else{
      if("nu" %in% names(params.list)){
        message("cor.fn = 'exponential'. Ignoring candidate values of nu.")
      }
      params.list[["nu"]] <- c(0.0)
    }

    if(!"boundary" %in% names(params.list)){
      message("Candidate values of boundary not specified. Using
              defaults c(0.5, 0.75).")
      params.list[["boundary"]] <- c(0.5, 0.75)
    }

    params.list <- params.list[c("phi", "nu", "boundary")]

  }

  # setup parameters for candidate models based on cartesian product of
  # candidate values of each parameter
  list_candidate <- candidate_models(params.list)

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
    }
    loopd.nMC <- loopd.controls[["nmc"]]
    if(loopd.nMC < 500){
      message("Number of Monte Carlo samples too low. Using defaults = 500.")
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
                    .Call("spGLMexactLOO", y, X, p, n, family, n.binom,
                          coords.D, cor.fn, V.beta, nu.beta, nu.z, sigmaSq.xi,
                          as.numeric(list_candidate[[x]]["phi"]),
                          as.numeric(list_candidate[[x]]["nu"]),
                          as.numeric(list_candidate[[x]]["boundary"]),
                          n.samples, loopd, loopd.method, CV.K, loopd.nMC,
                          verbose_child)
                          }, future.seed = TRUE)

  }else{

    # Get current plan invoked by future::plan() by the user
    current_plan <- future::plan()
    if(!inherits(current_plan, "sequential")){
      message("Parallelization plan other than 'sequential' setup but parallel
      is set to FALSE. Ignoring parallelization plan.")
    }

    samps <- lapply(1:length(list_candidate), function(x){
                    .Call("spGLMexactLOO", y, X, p, n, family, n.binom,
                          coords.D, cor.fn, V.beta, nu.beta, nu.z, sigmaSq.xi,
                          as.numeric(list_candidate[[x]]["phi"]),
                          as.numeric(list_candidate[[x]]["nu"]),
                          as.numeric(list_candidate[[x]]["boundary"]),
                          n.samples, loopd, loopd.method, CV.K, loopd.nMC,
                          verbose_child)
                        })

  }

  loopd_mat <- do.call("cbind", lapply(samps, function(x) x[["loopd"]]))
  # assists numerical stability for CVXR objective function evaluation
  # loopd_mat[loopd_mat < -10] <- -10
  # return(loopd_mat)

  out_CVXR <- get_stacking_weights(loopd_mat, solver = solver)
  # out_CVXR <- stacking_weights(loopd_mat, solver = solver)
  # w_hat <- loo::stacking_weights(loopd_mat)
  run.time <- proc.time() - ptm

  w_hat <- out_CVXR$weights
  w_hat <- as.numeric(w_hat)
  solver_status <- out_CVXR$status
  w_hat <- sapply(w_hat, function(x) max(0, x))
  w_hat <- w_hat / sum(w_hat)
  # solver_status <- "BFGS"

  stack_out <- as.matrix(do.call("rbind", lapply(list_candidate, unlist)))
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
  out$candidate.models <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$solver.status <- solver_status

  class(out) <- "spGLMstack"

  return(out)

}