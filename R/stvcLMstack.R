#' Bayesian spatially-temporally varying coefficients linear model using
#' predictive stacking
#'
#' @description Fits Bayesian linear models with spatially-temporally varying
#' coefficients for a Gaussian response on a collection of candidate models,
#' constructed from candidate values of the spatial-temporal process parameters
#' and the noise-to-spatial variance ratios supplied by the user, and combines
#' inference by stacking predictive densities.
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially-temporally varying
#' coefficients. See examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{stvcLMstack} is called.
#' @param sp_coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' spatial coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param time_coords an \eqn{n \times 1}{n x 1} matrix of the observation
#' temporal coordinates in \eqn{\mathcal{T} \subseteq [0, \infty)}.
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial-temporal dependence structure among the observations.
#' Supported covariance model key words are: \code{'gneiting-decay'} (Gneiting
#' and Guttorp 2010).
#' @param process.type a quoted keyword specifying the model for the
#' spatial-temporal processes of the varying coefficients: `'independent'` or
#' `'independent.shared'`. See [stvcLMexact()].
#' @param priors either \code{"flat"} (default), which assigns the prior
#' \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}, or a list with tags
#' \code{beta.norm} (a list containing \eqn{\mu_\beta} and \eqn{V_\beta})
#' and/or \code{sigma.sq.ig} (a vector containing \eqn{a_\sigma} and
#' \eqn{b_\sigma}). A component not supplied in the list receives its flat
#' prior.
#' @param candidate.models an object of class `candidateModels` containing a
#' list of candidate models for stacking, each with tags `phi_s`, `phi_t` and
#' `noise_sp_ratio`: vectors of length \eqn{r} if `process.type =
#' 'independent'` (use `list()` entries in [candidateModels()]), otherwise
#' scalars. See [candidateModels()] for details.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd.method character. Valid inputs are `'exact'` (default) and
#' `'PSIS'`. The option `'exact'` finds the exact leave-one-out predictive
#' densities in closed form. The option `'PSIS'` uses Pareto-smoothed
#' importance sampling (Vehtari *et al.* 2024); with many latent effects its
#' Pareto \eqn{k} diagnostics are often high, so `'exact'` is recommended.
#' @param parallel logical. If \code{parallel=FALSE}, the parallelization plan,
#' if set up by the user, is ignored. If \code{parallel=TRUE}, the function
#' inherits the parallelization plan that is set by the user via the function
#' [future::plan()] only. Depending on the parallel backend available, users
#' may choose their own plan. More details are available at
#' \url{https://cran.R-project.org/package=future}.
#' @param solver (optional) Specifies the name of the solver that will be used
#' to obtain optimal stacking weights for each candidate model. Default order
#' is \code{c("CLARABEL", "ECOS", "SCS")}. Users can use other solvers
#' supported by the \link[CVXR]{CVXR-package} package.
#' @param verbose logical. If \code{TRUE}, prints model-specific optimal
#' stacking weights.
#' @param ... currently no additional argument.
#' @return An object of class \code{stvcLMstack}, which is a list including the
#' following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 4, containing posterior
#'  samples of fixed effects (\code{beta}), the noise variance
#'  (\code{sigmaSq}), the process variances (\code{sigmaSq.z}) and the
#'  spatial-temporal effects (\code{z}) for that model, as in
#'  [stvcLMexact()].}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`model.params`}{a list with one element per candidate model, each a
#'  named list of its parameters: \code{phi_s}, \code{phi_t} and
#'  \code{noise_sp_ratio} (vectors of length \eqn{r} if \code{process.type =
#'  'independent'}).}
#' \item{`stacking.summary`}{a matrix with one row per candidate model,
#'  containing its parameters and its optimal stacking weight, for display.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#' candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`diagnostics`}{a list of diagnostics. Element \code{numerical} is a
#' data frame with one row per candidate model (per candidate model and
#' process, with columns \code{model} and \code{process}, if
#' \code{process.type = 'independent'}) and columns \code{min.pivot},
#' \code{min.cor} and \code{max.cor}, as in [stvcLMexact()]. If
#' \code{loopd.method = 'PSIS'}, element \code{pareto} contains the Pareto
#' \eqn{k} diagnostic values of each candidate model (\code{k}), the threshold
#' above which they are unreliable (\code{threshold}) and the number of values
#' above it for each model (\code{n.high}). Element \code{solver} describes the
#' optimization for the stacking weights: the solver used (\code{used}) and
#' its status (\code{status}), the installed and requested solvers, the search
#' order, the attempts with their status, and whether the fallback
#' \code{loo::stacking_weights()} was used. If \code{verbose = TRUE}, a
#' "Diagnostics" section is printed if there is an issue: numerical flags of
#' candidate models with stacking weight above 0.05 (extreme candidates with
#' negligible weight are expected in a stacking grid and are only counted),
#' Pareto \eqn{k} values above the threshold for any model, and solver
#' problems.}
#' }
#' This object can be used to make predictions at new locations or times with
#' [posteriorPredict()] and to sample from the stacked posterior with
#' [stackedSampler()].
#' @details Instead of assigning priors on the process parameters
#' \eqn{\phi_s}, \eqn{\phi_t} and the noise-to-spatial variance ratios
#' \eqn{\delta^2}, we consider a set of candidate models
#' \eqn{\mathcal{M} = \{M_1, \ldots, M_G\}} based on candidate values of these
#' parameters. For each \eqn{g}, we sample exactly from the posterior
#' distribution \eqn{p(\sigma^2, \beta, z \mid y, M_g)} (see [stvcLMexact()])
#' and find the leave-one-out predictive densities \eqn{p(y_i \mid y_{-i},
#' M_g)}. The stacking weights solve
#'  \deqn{
#'  \begin{aligned}
#'  \max_{w_1, \ldots, w_G}& \, \frac{1}{n} \sum_{i = 1}^n \log \sum_{g = 1}^G
#'  w_g p(y_i \mid y_{-i}, M_g) \\
#'  \text{subject to} & \quad w_g \geq 0, \sum_{g = 1}^G w_g = 1.
#'  \end{aligned}
#'  }
#' Candidate models that share the process parameters \eqn{(\phi_s, \phi_t)}
#' are fitted together, building and factorizing the spatial-temporal
#' correlation matrices only once.
#' @seealso [stvcLMexact()], [stvcGLMstack()], [spLMstack()]
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
#' \donttest{
#' set.seed(1234)
#' n <- 100
#' dat <- data.frame(s1 = runif(n), s2 = runif(n), t_coords = runif(n),
#'                   x1 = rnorm(n))
#' dat$slope <- 0.5 + sin(2 * pi * dat$s1) * cos(pi * dat$t_coords)
#' dat$y <- 1 + dat$slope * dat$x1 + rnorm(n, sd = 0.5)
#'
#' # processes sharing (phi_s, phi_t, noise_sp_ratio): scalar candidates
#' mod.list <- candidateModels(list(phi_s = c(1, 3), phi_t = c(0.5, 2),
#'                                  noise_sp_ratio = c(0.5, 2)), "cartesian")
#' mod1 <- stvcLMstack(y ~ x1 + (x1), data = dat,
#'                     sp_coords = as.matrix(dat[, c("s1", "s2")]),
#'                     time_coords = as.matrix(dat[, "t_coords"]),
#'                     cor.fn = "gneiting-decay",
#'                     process.type = "independent.shared",
#'                     candidate.models = mod.list,
#'                     n.samples = 500)
#' post_samps <- stackedSampler(mod1)
#' slope <- sweep(post_samps$z[n + 1:n, ], 2, post_samps$beta[2, ], "+")
#' cor(apply(slope, 1, median), dat$slope)
#'
#' # independent processes (r = 2): vector-valued candidates via list()
#' mod.list2 <- candidateModels(list(phi_s = list(c(1, 1), c(3, 3)),
#'                                   phi_t = list(c(1, 1)),
#'                                   noise_sp_ratio = list(c(1, 1), c(2, 0.5))),
#'                              "cartesian")
#' mod2 <- stvcLMstack(y ~ x1 + (x1), data = dat,
#'                     sp_coords = as.matrix(dat[, c("s1", "s2")]),
#'                     time_coords = as.matrix(dat[, "t_coords"]),
#'                     cor.fn = "gneiting-decay",
#'                     process.type = "independent",
#'                     candidate.models = mod.list2,
#'                     n.samples = 500, verbose = FALSE)
#' }
#' @export
stvcLMstack <- function(formula, data = parent.frame(), sp_coords, time_coords,
                        cor.fn, process.type, priors = "flat",
                        candidate.models, n.samples, loopd.method = "exact",
                        parallel = FALSE, solver = NULL, verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

  ##### process type #####
  process.type <- check_stvc_lm_process_type(process.type)

  ##### formula and coordinates #####
  if(missing(sp_coords)){
    stop("sp_coords must be supplied.")
  }
  if(missing(time_coords)){
    stop("time_coords must be supplied.")
  }
  dd <- stvc_lm_data(formula, data, sp_coords, time_coords)
  n <- dd$n
  p <- dd$p
  r <- dd$r
  nR <- if(process.type == "independent") r else 1L

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("cor.fn must be specified")
  }
  if(!identical(cor.fn, "gneiting-decay")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from c('gneiting-decay')")
  }

  ##### priors #####
  pr <- parse_lm_priors(priors, p)

  ## sample size check: a flat prior on beta requires n - 1 > p for the
  ## leave-one-out predictive densities
  if(pr$beta.prior == "flat" && n <= p + 1){
    stop("a flat prior on beta requires n - 1 > p for leave-one-out predictive densities; supply beta.norm in priors.")
  }

  #### candidate models ####
  if(missing(candidate.models)){
    stop("candidate.models must be supplied.")
  }
  if(!inherits(candidate.models, "candidateModels")){
    stop("candidate.models must be an object of class 'candidateModels'.")
  }
  list_candidate <- candidate.models
  check_validity <- all(vapply(list_candidate, function(x){
    is.list(x) && length(x) == 3 &&
      identical(sort(names(x)), c("noise_sp_ratio", "phi_s", "phi_t")) &&
      all(vapply(x, function(v) is.numeric(v) && length(v) == nR && all(is.finite(v)) && all(v > 0),
                 logical(1)))
  }, logical(1)))
  if(!check_validity){
    if(nR > 1){
      stop("each element of candidate.models must be a named list with entries 'phi_s', 'phi_t' and 'noise_sp_ratio', each a vector of ", nR,
           " positive reals when process.type = 'independent' (supply them as list() entries in candidateModels()).")
    }
    stop("each element of candidate.models must be a named list with positive scalar entries 'phi_s', 'phi_t' and 'noise_sp_ratio'.")
  }

  #### Leave-one-out setup ####
  loopd <- 1L
  loopd.method <- tolower(loopd.method)
  if(!loopd.method %in% c("exact", "psis")){
    stop("loopd.method = '", loopd.method, "' is not a valid option; choose from c('exact', 'PSIS').")
  }

  ##### sampling setup #####
  if(missing(n.samples)){
    stop("n.samples must be specified.")
  }
  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  # candidate models sharing (phi_s, phi_t) are fitted in one call, which builds and factorizes the
  # spatial-temporal correlation matrices once and loops over their noise_sp_ratio values
  cand_phi_s <- lapply(list_candidate, function(x) as.double(x[["phi_s"]]))
  cand_phi_t <- lapply(list_candidate, function(x) as.double(x[["phi_t"]]))
  cand_deltasq <- lapply(list_candidate, function(x) as.double(x[["noise_sp_ratio"]]))
  cand_key <- vapply(seq_along(list_candidate), function(x){
    paste(sprintf("%a", c(cand_phi_s[[x]], cand_phi_t[[x]])), collapse = " ")   # exact (hexadecimal) keys
  }, character(1))
  cand_groups <- split(seq_along(list_candidate), factor(cand_key, levels = unique(cand_key)))
  names(cand_groups) <- NULL

  # fits the candidate models in group g
  fit_group <- function(g){
    idx <- cand_groups[[g]]
    .Call(C_stvcLMexactGrid, dd$y, dd$X, dd$X_tilde, n, p, r,
          dd$sp_coords, dd$time_coords, cor.fn, process.type,
          cand_phi_s[[idx[1]]], cand_phi_t[[idx[1]]],
          pr$beta.prior, pr$beta.Norm, pr$sigma.sq.IG,
          unlist(cand_deltasq[idx]), n.samples, loopd, loopd.method)
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
                 NWORKERS.machine - 1, "to avoid overcommitment of resources."))
    }

    if(rstudioapi::isAvailable()){
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

    samps <- ungroup(future_lapply(seq_along(cand_groups), fit_group, future.seed = TRUE))

  }else{

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
  out <- suppressMessages(get_stacking_weights(loopd_mat, solver = solver, verbose = FALSE))

  w_hat <- out$weights
  if(identical(out$solver, "none")){
    message(loo_fallback_message())
  }
  solver_status <- out$status
  solver_used <- out$solver
  out_solver_details <- out$details

  run.time <- proc.time() - ptm

  # columns in a fixed order: phi_s, phi_t, noise_sp_ratio (one column per process if independent)
  stack_out <- as.matrix(do.call("rbind", lapply(list_candidate, function(x){
    unlist(lapply(c("phi_s", "phi_t", "noise_sp_ratio"), function(nm) as.numeric(x[[nm]])))
  })))
  par_names <- if(nR > 1){
    unlist(lapply(c("phi_s", "phi_t", "noise_sp_ratio"), function(nm) paste0(nm, "[", seq_len(nR), "]")))
  }else{
    c("phi_s", "phi_t", "noise_sp_ratio")
  }
  stack_out <- cbind(stack_out, round(w_hat, 3))
  colnames(stack_out) <- c(par_names, "weight")
  rownames(stack_out) <- paste("Model", seq_len(nrow(stack_out)))

  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", seq_along(list_candidate), sep = "")

  diagnostics <- list(numerical = collect_diagnostics(samps, paste("Model", seq_along(list_candidate))))
  if(loopd.method == "psis"){
    pareto_k_list <- lapply(samps, function(x) x[["loopd.pareto_k"]])
    names(pareto_k_list) <- names(loopd_list)
    diagnostics$pareto <- pareto_diagnostics(pareto_k_list, n.samples)
  }
  diagnostics$solver <- c(list(used = solver_used, status = solver_status), out_solver_details)
  if(verbose){
    print_diagnostics(diagnostics, weights = w_hat,
                      pivot.hint = "nearly coincident space-time locations, very small decay parameters phi_s, phi_t, or very small noise-to-spatial variance ratios")
  }

  samps <- lapply(samps, function(x) x[c("beta", "sigmaSq", "sigmaSq.z", "z")])
  names(samps) <- paste("Model", seq_along(list_candidate), sep = "")

  out <- list()
  out$y <- dd$y
  out$X <- dd$X
  out$X.names <- dd$X.names
  out$X.stvc.names <- dd$X_tilde.names
  out$sp_coords <- dd$sp_coords
  out$time_coords <- dd$time_coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$priors <- pr$out
  out$n.samples <- n.samples
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.method
  out$n.models <- length(list_candidate)
  # model parameters of each candidate, read by name (e.g. by posteriorPredict())
  out$model.params <- lapply(seq_along(list_candidate), function(k){
    list(phi_s = cand_phi_s[[k]], phi_t = cand_phi_t[[k]], noise_sp_ratio = cand_deltasq[[k]])
  })
  names(out$model.params) <- paste("Model", seq_along(list_candidate))
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$diagnostics <- diagnostics

  class(out) <- "stvcLMstack"

  return(out)

}
