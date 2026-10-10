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
#' @param priors either \code{"flat"} (default), which assigns the prior
#'  \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}, or a list with tags
#'  \code{beta.norm} (a list containing \eqn{\mu_\beta} and \eqn{V_\beta})
#'  and/or \code{sigma.sq.ig} (a vector containing \eqn{a_\sigma} and
#'  \eqn{b_\sigma}). A component not supplied in the list receives its flat
#'  prior, \eqn{p(\beta) \propto 1} or \eqn{p(\sigma^2) \propto 1/\sigma^2}.
#' @param candidate.models an object of class `candidateModels` containing a
#'  list of candidate models for stacking. See [candidateModels()] for details.
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
#'  to obtain optimal stacking weights for each candidate model. Default order
#'  is \code{c("CLARABEL", "ECOS", "SCS")}. Users can use other solvers
#'  supported by the \link[CVXR]{CVXR-package} package.
#' @param verbose logical. If \code{TRUE}, prints model-specific optimal
#'  stacking weights.
#' @param ... currently no additional argument.
#' @return An object of class \code{spLMstack}, which is a list including the
#'  following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 4, containing posterior
#'  samples of fixed effects (\code{beta}), measurement error variance
#'  (\code{sigmaSq}), spatial variance (\code{sigmaSq.z}), and spatial effects
#'  (\code{z}) for that model.}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`model.params`}{a list with one element per candidate model, each a
#'  named list of its parameters: \code{phi}, \code{nu} (\code{NA} for the
#'  exponential correlation function) and \code{noise_sp_ratio} (noise-to-spatial variance ratio).}
#' \item{`stacking.summary`}{a matrix with one row per candidate model,
#'  containing its parameters and its optimal stacking weight, for display.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#'  candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`diagnostics`}{a list of diagnostics. Element \code{numerical} is a
#'  data frame with one row per candidate model and columns \code{min.pivot}
#'  (the smallest relative Cholesky pivot of the \eqn{n \times n}{n x n}
#'  factorizations; values below 1e-8 indicate a nearly singular covariance
#'  matrix), \code{min.cor} and \code{max.cor} (the correlations of the two
#'  farthest-apart and of the two closest locations; values of \code{min.cor}
#'  above 0.95 suggest an effective range far exceeding the extent of the data,
#'  values of \code{max.cor} below 0.05 nearly uncorrelated locations),
#'  obtained from quantities the fits compute anyway. If
#'  \code{loopd.method = 'PSIS'}, element \code{pareto} is a list with the
#'  Pareto \eqn{k} diagnostic values of each candidate model (\code{k}), the
#'  threshold above which they are unreliable (\code{threshold}) and the
#'  number of values above it for each model (\code{n.high}). Element
#'  \code{solver} describes the optimization for the stacking weights: the
#'  solver used (\code{used}) and its status (\code{status}), the installed
#'  and requested solvers, the search order, the attempts with their status,
#'  and whether the fallback \code{loo::stacking_weights()} was used. If
#'  \code{verbose = TRUE}, a "Diagnostics" section is printed if there is an
#'  issue: numerical flags of candidate models with stacking weight above 0.05
#'  (extreme candidates with negligible weight are expected in a stacking grid
#'  and are only counted), Pareto \eqn{k} values above the threshold for any
#'  model, and solver problems (a requested solver not installed, an
#'  inaccurate solution, or the fallback).}
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
#' cand.mod <- candidateModels(list(phi = c(1.5, 3),
#'                                  nu = c(0.5, 1),
#'                                  noise_sp_ratio = c(1)),
#'                             "cartesian")
#'
#' mod1 <- spLMstack(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   priors = prior_list,
#'                   candidate.models = cand.mod,
#'                   n.samples = 1000, loopd.method = "exact",
#'                   parallel = FALSE, verbose = TRUE)
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
                      priors = "flat", candidate.models, n.samples,
                      loopd.method,
                      parallel = FALSE, solver = NULL, verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

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

  check_no_missing(y = y, X = X, coords = coords)
  check_distinct_coords(coords)

  ## distances are computed in C++ from the coordinates
  storage.mode(coords) <- "double"

  ##### correlation function #####
  if (missing(cor.fn)) {
    stop("error: cor.fn must be specified")
  }
  if (!cor.fn %in% c("exponential", "matern")) {
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
    c('exponential', 'matern').")
  }

  ##### priors #####
  # priors = "flat" assigns p(beta, sigma.sq) proportional to 1/sigma.sq; if
  # priors is a list, a component not supplied receives its flat prior, i.e.,
  # p(beta) proportional to 1 or p(sigma.sq) proportional to 1/sigma.sq
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.prior <- "flat"
  sigma.sq.IG <- c(0.0, 0.0)

  if(is.character(priors)){

    if(length(priors) != 1 || tolower(priors) != "flat"){
      stop("error: priors must be either 'flat' or a named list with tags
           'beta.norm' and/or 'sigma.sq.ig'.")
    }

  }else if(is.list(priors)){

    if(is.null(names(priors))){
      stop("error: priors must be either 'flat' or a named list with tags
           'beta.norm' and/or 'sigma.sq.ig'.")
    }
    names(priors) <- tolower(names(priors))
    if(any(!names(priors) %in% c("beta.norm", "sigma.sq.ig"))){
      stop("error: invalid tag(s) in priors: '",
           paste(setdiff(names(priors), c("beta.norm", "sigma.sq.ig")),
                 collapse = "', '"),
           "'. Valid tags are 'beta.norm' and 'sigma.sq.ig'.")
    }

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
                   " covariance matrix.", sep = ""))
      }
      check_cov_matrix(beta.Norm[[2]], p, "the prior covariance of beta (beta.norm[[2]])")
      storage.mode(beta.Norm[[1]]) <- "double"
      storage.mode(beta.Norm[[2]]) <- "double"
      beta.prior <- "normal"
    }

    ## Setup prior for sigma.sq
    if("sigma.sq.ig" %in% names(priors)){
      sigma.sq.IG <- priors[["sigma.sq.ig"]]
      if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){
        stop("error: sigma.sq.IG must be a vector of length 2")
      }
      if(any(sigma.sq.IG <= 0)){
        stop("error: sigma.sq.IG must be a positive vector of length 2")
      }
      sigma.sq.prior <- "ig"
    }

  }else{
    stop("error: priors must be either 'flat' or a named list with tags
         'beta.norm' and/or 'sigma.sq.ig'.")
  }

  ## storage mode
  storage.mode(sigma.sq.IG) <- "double"

  ## sample size check: a flat prior on beta requires n - 1 > p for the
  ## leave-one-out predictive densities
  if(beta.prior == "flat"){
    if(n <= p + 1){
      stop("error: a flat prior on beta requires n - 1 > p for leave-one-out
           predictive densities; supply beta.norm in priors.")
    }
  }

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
        identical(sort(names(x)), c("noise_sp_ratio", "nu", "phi")) &&
        all(vapply(x, function(v) is.numeric(v) && length(v) == 1, logical(1)))
      }, logical(1)))
      if(!check_validity){
        stop("error: each element of candidate.models must be a named list with
             scalar numeric entries 'phi', 'nu' and 'noise_sp_ratio'.")
      }
    }else{
      check_validity <- all(vapply(candidate.models, function(x){
        length(x) == 2 &&
        identical(sort(names(x)), c("noise_sp_ratio", "phi")) &&
        all(vapply(x, function(v) is.numeric(v) && length(v) == 1, logical(1)))
      }, logical(1)))
      if(!check_validity){
        stop("error: each element of candidate.models must be a named list with
             scalar numeric entries 'phi' and 'noise_sp_ratio'.")
      }
      candidate.models <- lapply(candidate.models, function(x){
        x[["nu"]] <- NA_real_                                  # not used by the exponential
        x
      })
      class(candidate.models) <- "candidateModels"
    }
  }

  list_candidate <- candidate.models

  #### Leave-one-out setup ####
  loopd <- TRUE

  if(missing(loopd.method)){
    message("loopd.method not specified. Using 'exact'.")
    loopd.method <- "exact"
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

  # candidate models sharing (phi, nu) are fitted in one call, which builds the
  # spatial correlation matrix once and loops over their noise_sp_ratio values
  cand_phi <- vapply(list_candidate, function(x) as.numeric(x[["phi"]]), numeric(1))
  cand_nu <- vapply(list_candidate, function(x) as.numeric(x[["nu"]]), numeric(1))
  cand_deltasq <- vapply(list_candidate, function(x) as.numeric(x[["noise_sp_ratio"]]),
                         numeric(1))
  cand_key <- sprintf("%a %a", cand_phi, cand_nu)          # exact (hexadecimal) keys
  cand_groups <- split(seq_along(list_candidate),
                       factor(cand_key, levels = unique(cand_key)))
  names(cand_groups) <- NULL

  # fits the candidate models in group g
  fit_group <- function(g){
    idx <- cand_groups[[g]]
    .Call(C_spLMexactLOOgrid, y, X, p, n, coords,
          beta.prior, beta.Norm, sigma.sq.IG,
          cand_phi[idx[1]], cand_nu[idx[1]], cand_deltasq[idx],
          cor.fn, n.samples, loopd, loopd.method)
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
    unlist(x[c("phi", "nu", "noise_sp_ratio")])
  })))
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

  diagnostics <- list(numerical = collect_diagnostics(samps, paste("Model", seq_along(list_candidate))))
  if(loopd.method == "psis"){
    pareto_k_list <- lapply(samps, function(x) x[["loopd.pareto_k"]])
    names(pareto_k_list) <- names(loopd_list)
    diagnostics$pareto <- pareto_diagnostics(pareto_k_list, n.samples)
  }
  diagnostics$solver <- c(list(used = solver_used, status = solver_status), out_solver_details)
  if(verbose){
    print_diagnostics(diagnostics, weights = w_hat,
                      pivot.hint = "nearly coincident locations, a very small decay parameter, or a very small noise-to-spatial variance ratio")
  }

  samps <- lapply(samps, function(x) x[c("beta", "sigmaSq", "sigmaSq.z", "z")])
  names(samps) <- paste("Model", 1:length(list_candidate), sep = "")

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  if(beta.prior == "normal"){
    beta.Norm.out <- list(mu = beta.Norm[[1]], V = matrix(beta.Norm[[2]], p, p))
  }else{
    beta.Norm.out <- "flat"
  }
  if(sigma.sq.prior == "ig"){
    sigma.sq.IG.out <- sigma.sq.IG
  }else{
    sigma.sq.IG.out <- "flat"
  }
  out$priors <- list(beta.Norm = beta.Norm.out, sigma.sq.IG = sigma.sq.IG.out)
  out$n.samples <- n.samples
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.method
  out$n.models <- length(list_candidate)
  # model parameters of each candidate, read by name (e.g. by posteriorPredict())
  out$model.params <- lapply(list_candidate, function(x){
    list(phi = as.numeric(x[["phi"]]), nu = as.numeric(x[["nu"]]),
         noise_sp_ratio = as.numeric(x[["noise_sp_ratio"]]))
  })
  names(out$model.params) <- paste("Model", seq_along(list_candidate))
  # table of the candidate models and their stacking weights (for display)
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$diagnostics <- diagnostics

  class(out) <- "spLMstack"

  return(out)

}