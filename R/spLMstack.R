#' Bayesian spatial linear model using predictive stacking
#'
#' @description Fits Bayesian spatial linear model on a collection of candidate
#' models constructed based on some candidate values of some model parameters
#' specified by the user and subsequently combines inference by stacking
#' predictive densities. See Zhang, Tang and Banerjee (2025) for more details.
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#'  If not found in \code{data}, the variables are taken from
#'  \code{environment(formula)}, typically the environment from which
#'  \code{spLMstack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors a list with each tag corresponding to a parameter name and
#'  containing prior details. If not supplied, uses defaults.
#' @param params.list a list containing candidate values of spatial process
#'  parameters for the `cor.fn` used, and, noise-to-spatial variance ratio.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd.method character. Valid inputs are `'exact'` and `'PSIS'`. The
#'  option `'exact'` corresponds to exact leave-one-out predictive densities.
#'  The option `'PSIS'` is faster, as it finds approximate leave-one-out
#'  predictive densities using Pareto-smoothed importance sampling
#'  (Gelman *et al.* 2024).
#' @param parallel logical. If \code{parallel=FALSE}, the parallelization plan,
#'  if set up by the user, is ignored. If \code{parallel=TRUE}, the function
#'  inherits the parallelization plan that is set by the user via the function
#'  [future::plan()] only. Depending on the parallel backend available, users
#'  may choose their own plan. More details are available at
#'  \url{https://cran.R-project.org/package=future}.
#' @param solver (optional) Specifies the name of the solver that will be used
#'  to obtain optimal stacking weights for each candidate model. Default is
#'  \code{"ECOS"}. Users can use other solvers supported by the
#'  \link[CVXR]{CVXR-package} package.
#' @param verbose logical. If \code{TRUE}, prints model-specific optimal
#'  stacking weights.
#' @param ... currently no additional argument.
#' @return An object of class \code{spLMstack}, which is a list including the
#'  following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 3, containing posterior
#'  samples of fixed effects (\code{beta}), variance parameter
#'  (\code{sigmaSq}), spatial effects (\code{z}) for that model.}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
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
#'  and \eqn{\nu}, noise-to-spatial variance ratio \eqn{\delta^2}, we consider
#'  a set of candidate models based on some candidate values of these parameters
#'  supplied by the user. Suppose the set of candidate models is
#'  \eqn{\mathcal{M} = \{M_1, \ldots, M_G\}}. Then for each
#'  \eqn{g = 1, \ldots, G}, we sample from the posterior distribution
#'  \eqn{p(\sigma^2, \beta, z \mid y, M_g)} under the model \eqn{M_g} and find
#'  leave-one-out predictive densities \eqn{p(y_i \mid y_{-i}, M_g)}. Then we
#'  solve the optimization problem
#'  \deqn{
#'  \begin{aligned}
#'  \max_{w_1, \ldots, w_G}& \, \frac{1}{n} \sum_{i = 1}^n \log \sum_{g = 1}^G
#'  w_g p(y_i \mid y_{-i}, M_g) \\
#'  \text{subject to} & \quad w_g \geq 0, \sum_{g = 1}^G w_g = 1
#'  \end{aligned}
#'  }
#' to find the optimal stacking weights \eqn{\hat{w}_1, \ldots, \hat{w}_G}.
#' @seealso [spLMexact()], [spGLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto
#'  Smoothed Importance Sampling." *Journal of Machine Learning Research*,
#'  **25**(72), 1-58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @references Zhang L, Tang W, Banerjee S (2025). "Bayesian Geostatistics Using
#' Predictive Stacking." *Journal of the American Statistical Association*,
#' **In press**. \doi{10.1080/01621459.2025.2566449}.
#' @importFrom rstudioapi isAvailable
#' @importFrom parallel detectCores
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_lapply
#' @examples
#' set.seed(1234)
#' # load data and work with first 100 rows
#' data(simGaussian)
#' dat <- simGaussian[1:100, ]
#'
#' # setup prior list
#' muBeta <- c(0, 0)
#' VBeta <- cbind(c(1.0, 0.0), c(0.0, 1.0))
#' sigmaSqIGa <- 2
#' sigmaSqIGb <- 2
#' prior_list <- list(beta.norm = list(muBeta, VBeta),
#'                    sigma.sq.ig = c(sigmaSqIGa, sigmaSqIGb))
#'
#' mod1 <- spLMstack(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   priors = prior_list,
#'                   params.list = list(phi = c(1.5, 3),
#'                                      nu = c(0.5, 1),
#'                                      noise_sp_ratio = c(1)),
#'                   n.samples = 1000, loopd.method = "exact",
#'                   parallel = FALSE, solver = "ECOS", verbose = TRUE)
#'
#' post_samps <- stackedSampler(mod1)
#' post_beta <- post_samps$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' post_z <- post_samps$z
#' post_z_summ <- t(apply(post_z, 1,
#'                        function(x) quantile(x, c(0.025, 0.5, 0.975))))
#'
#' z_combn <- data.frame(z = dat$z_true,
#'                       zL = post_z_summ[, 1],
#'                       zM = post_z_summ[, 2],
#'                       zU = post_z_summ[, 3])
#'
#' library(ggplot2)
#' plot1 <- ggplot(data = z_combn, aes(x = z)) +
#'   geom_point(aes(y = zM), size = 0.25,
#'              color = "darkblue", alpha = 0.5) +
#'   geom_errorbar(aes(ymin = zL, ymax = zU),
#'                 width = 0.05, alpha = 0.15) +
#'   geom_abline(slope = 1, intercept = 0,
#'               color = "red", linetype = "solid") +
#'   xlab("True z") + ylab("Stacked posterior of z") +
#'   theme_bw() +
#'   theme(panel.background = element_blank(),
#'         aspect.ratio = 1)
#' @export
spLMstack <- function(formula, data = parent.frame(), coords, cor.fn,
                      priors, params.list, n.samples, loopd.method,
                      parallel = FALSE, solver = "ECOS", verbose = TRUE, ...){

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
    beta.prior <- "normal"
    beta.Norm <- list(rep(0.0, p), diag(100.0, p))
    sigma.sq.IG <- c(2, 2)
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
    if(!"sigma.sq.ig" %in% names(priors)){
      stop("error: sigma.sq.IG must be specified")
    }
    sigma.sq.IG <- priors[["sigma.sq.ig"]]

    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){
      stop("error: sigma.sq.IG must be a vector of length 2")
    }
    if(any(sigma.sq.IG <= 0)){
      stop("error: sigma.sq.IG must be a positive vector of length 2")
    }

  }

  ## storage mode
  storage.mode(sigma.sq.IG) <- "double"

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

    if(!"noise_sp_ratio" %in% names(params.list)){
      message("Candidate values of noise_sp_ratio not specified. Using
              defaults c(0.25, 1, 2).")
      params.list[["noise_sp_ratio"]] <- c(0.25, 1.0, 2.0)
    }

    params.list <- params.list[c("phi", "nu", "noise_sp_ratio")]

  }

  # setup parameters for candidate models based on cartesian product of
  # candidate values of each parameter
  list_candidate <- candidate_models(params.list)

  #### Leave-one-out setup ####
  loopd <- TRUE

  if(missing(loopd.method)){
    message("loopd.method not specified. Using 'exact'.")
  }

  loopd.method <- tolower(loopd.method)

  if(!loopd.method %in% c("exact", "psis")){
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
                          cor.fn, n.samples, loopd, loopd.method,
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
                    .Call("spLMexactLOO", y, X, p, n, coords.D,
                          beta.prior, beta.Norm, sigma.sq.IG,
                          as.numeric(list_candidate[[x]]["phi"]),
                          as.numeric(list_candidate[[x]]["nu"]),
                          as.numeric(list_candidate[[x]]["noise_sp_ratio"]),
                          cor.fn, n.samples, loopd, loopd.method,
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
  colnames(stack_out) = c("phi", "nu", "noise_sp_ratio", "weight")
  rownames(stack_out) = paste("Model", 1:nrow(stack_out))

  if(cor.fn == 'exponential'){
    stack_out <- stack_out[, c("phi", "noise_sp_ratio", "weight")]
  }

  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", 1:length(list_candidate), sep = "")

  samps <- lapply(samps, function(x) x[1:3])
  names(samps) <- paste("Model", 1:length(list_candidate), sep = "")

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$priors <- list(beta.Norm = list(mu = beta.Norm[[1]],
                                      V = matrix(beta.Norm[[2]], p, p)),
                     sigma.sq.IG = sigma.sq.IG)
  out$n.samples <- n.samples
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.method
  out$n.models <- length(list_candidate)
  out$candidate.models <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$solver.status <- solver_status

  class(out) <- "spLMstack"

  return(out)

}