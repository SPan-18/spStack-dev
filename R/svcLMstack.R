#' Bayesian spatially varying coefficients linear model using predictive
#' stacking
#'
#' @description Fits Bayesian linear models with spatially varying coefficients
#' for a Gaussian response on a collection of candidate models, constructed
#' from candidate values of the spatial process parameters and the
#' noise-to-spatial variance ratios supplied by the user, and combines
#' inference by stacking predictive densities.
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially varying coefficients. See
#' examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{svcLMstack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial dependence structure among the observations. Supported
#' covariance model key words are: \code{'exponential'} and \code{'matern'}.
#' @param process.type a quoted keyword specifying the model for the spatial
#' processes of the varying coefficients: `'independent'` or
#' `'independent.shared'`. See [svcLMexact()].
#' @param priors either \code{"flat"} (default), which assigns the prior
#' \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}, or a list with tags
#' \code{beta.norm} (a list containing \eqn{\mu_\beta} and \eqn{V_\beta})
#' and/or \code{sigma.sq.ig} (a vector containing \eqn{a_\sigma} and
#' \eqn{b_\sigma}). A component not supplied in the list receives its flat
#' prior.
#' @param candidate.models an object of class `candidateModels` containing a
#' list of candidate models for stacking, each with tags `phi`, `nu` (only if
#' `cor.fn = 'matern'`) and `noise_sp_ratio`: vectors of length \eqn{r} if
#' `process.type = 'independent'` (use `list()` entries in
#' [candidateModels()]), otherwise scalars. See [candidateModels()] for
#' details.
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
#' @return An object of class \code{svcLMstack}, which is a list including the
#' following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 4, containing posterior
#'  samples of fixed effects (\code{beta}), the noise variance
#'  (\code{sigmaSq}), the process variances (\code{sigmaSq.z}) and the
#'  spatially varying coefficients (\code{z}) for that model, as in
#'  [svcLMexact()].}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`model.params`}{a list with one element per candidate model, each a
#'  named list of its parameters: \code{phi}, \code{nu} (\code{NA} for the
#'  exponential correlation function) and \code{noise_sp_ratio} (vectors of
#'  length \eqn{r} if \code{process.type = 'independent'}).}
#' \item{`stacking.summary`}{a matrix with one row per candidate model,
#'  containing its parameters and its optimal stacking weight, for display.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#' candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`diagnostics`}{a list of diagnostics. Element \code{numerical} is a
#' data frame with one row per candidate model (per candidate model and
#' process, with columns \code{model} and \code{process}, if
#' \code{process.type = 'independent'}) and columns \code{min.pivot},
#' \code{min.cor} and \code{max.cor}, as in [svcLMexact()]. If
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
#' This object can be used to make predictions at new locations with
#' [posteriorPredict()] and to sample from the stacked posterior with
#' [stackedSampler()].
#' @details Instead of assigning priors on the process parameters \eqn{\phi},
#' \eqn{\nu} and the noise-to-spatial variance ratios \eqn{\delta^2}, we
#' consider a set of candidate models \eqn{\mathcal{M} = \{M_1, \ldots, M_G\}}
#' based on candidate values of these parameters. For each \eqn{g}, we sample
#' exactly from the posterior distribution \eqn{p(\sigma^2, \beta, z \mid y,
#' M_g)} (see [svcLMexact()]) and find the leave-one-out predictive densities
#' \eqn{p(y_i \mid y_{-i}, M_g)}. The stacking weights solve
#'  \deqn{
#'  \begin{aligned}
#'  \max_{w_1, \ldots, w_G}& \, \frac{1}{n} \sum_{i = 1}^n \log \sum_{g = 1}^G
#'  w_g p(y_i \mid y_{-i}, M_g) \\
#'  \text{subject to} & \quad w_g \geq 0, \sum_{g = 1}^G w_g = 1.
#'  \end{aligned}
#'  }
#' Candidate models that share the process parameters \eqn{(\phi, \nu)} are
#' fitted together, building and factorizing the spatial correlation matrices
#' only once, and the groups of such candidates are the units of parallel
#' work.
#' @seealso [svcLMexact()], [svcGLMstack()], [spLMstack()], [stvcLMstack()]
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
#' data(simSpatial)
#' dat <- simSpatial[1:200, ]
#'
#' # independent processes (r = 2): vector-valued candidates via list()
#' mod.list <- candidateModels(list(phi = list(c(4, 4), c(8, 4), c(8, 8)),
#'                                  nu = list(c(1.5, 1.5)),
#'                                  noise_sp_ratio = list(c(0.1, 0.1), c(0.25, 0.25))),
#'                             "cartesian")
#' mod1 <- svcLMstack(y_svc ~ x1 + x2 + (x1), data = dat,
#'                    coords = as.matrix(dat[, c("s1", "s2")]),
#'                    cor.fn = "matern", process.type = "independent",
#'                    candidate.models = mod.list, n.samples = 500)
#' post_samps <- stackedSampler(mod1)
#' n <- nrow(dat)
#' cor(apply(post_samps$z[n + 1:n, ], 1, median), dat$z_true)
#' @export
svcLMstack <- function(formula, data = parent.frame(), coords, cor.fn,
                       process.type, priors = "flat", candidate.models,
                       n.samples, loopd.method = "exact", parallel = FALSE,
                       solver = NULL, verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

  ##### process type #####
  process.type <- check_stvc_lm_process_type(process.type)

  ##### formula and coordinates #####
  if(missing(coords)){
    stop("coords must be supplied.")
  }
  dd <- stvc_lm_data(formula, data, coords, NULL, coords.name = "coords")
  n <- dd$n
  p <- dd$p
  r <- dd$r
  nR <- vc_nR(process.type, r)

  ##### correlation function #####
  cor.fn <- check_svc_cor_fn(if(missing(cor.fn)) NULL else cor.fn)

  ##### priors #####
  pr <- parse_lm_priors(priors, p)

  ## sample size check: a flat prior on beta requires n - 1 > p for the
  ## leave-one-out predictive densities
  if(pr$beta.prior == "flat" && n <= p + 1){
    stop("a flat prior on beta requires n - 1 > p for leave-one-out predictive densities; supply beta.norm in priors.")
  }

  #### candidate models ####
  cands <- svc_candidates(if(missing(candidate.models)) NULL else candidate.models,
                          cor.fn, nR, process.type, "noise_sp_ratio")

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

  # candidate models sharing (phi, nu) are fitted in one call, which builds and
  # factorizes the spatial correlation matrices once and loops over their
  # noise_sp_ratio values
  cand_phi <- lapply(cands, function(x) x[["phi"]])
  cand_nu <- lapply(cands, function(x) svc_c_nu(x[["nu"]]))
  cand_deltasq <- lapply(cands, function(x) x[["noise_sp_ratio"]])
  cand_groups <- vc_candidate_groups(cand_phi, cand_nu)

  # fits the candidate models in group g (the routine of stvcLMstack() with a
  # spatial correlation function: no temporal coordinates, and nu in place of
  # phi_t)
  fit_group <- function(g){
    idx <- cand_groups[[g]]
    .Call(C_stvcLMexactGrid, dd$y, dd$X, dd$X_tilde, n, p, r,
          dd$sp_coords, NULL, cor.fn, process.type,
          cand_phi[[idx[1]]], cand_nu[[idx[1]]],
          pr$beta.prior, pr$beta.Norm, pr$sigma.sq.IG,
          unlist(cand_deltasq[idx]), n.samples, loopd, loopd.method)
  }

  #### main function call ####
  ptm <- proc.time()

  samps <- vc_run_groups(cand_groups, fit_group, length(cands), parallel)

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

  stack_out <- svc_stacking_summary(cands, cor.fn, nR, "noise_sp_ratio", w_hat)
  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  model_names <- paste("Model", seq_along(cands))
  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", seq_along(cands), sep = "")

  diagnostics <- list(numerical = collect_diagnostics(samps, model_names))
  if(loopd.method == "psis"){
    pareto_k_list <- lapply(samps, function(x) x[["loopd.pareto_k"]])
    names(pareto_k_list) <- names(loopd_list)
    diagnostics$pareto <- pareto_diagnostics(pareto_k_list, n.samples)
  }
  diagnostics$solver <- c(list(used = solver_used, status = solver_status), out_solver_details)
  if(verbose){
    print_diagnostics(diagnostics, weights = w_hat,
                      pivot.hint = "nearly coincident locations, very small decay parameters phi, or very small noise-to-spatial variance ratios")
  }

  samps <- lapply(samps, function(x) x[c("beta", "sigmaSq", "sigmaSq.z", "z")])
  names(samps) <- names(loopd_list)

  out <- list()
  out$y <- dd$y
  out$X <- dd$X
  out$X.names <- dd$X.names
  out$X.svc.names <- dd$X_tilde.names
  out$coords <- dd$sp_coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$priors <- pr$out
  out$n.samples <- n.samples
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- loopd.method
  out$n.models <- length(cands)
  # model parameters of each candidate, read by name (e.g. by posteriorPredict())
  out$model.params <- cands
  names(out$model.params) <- model_names
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$diagnostics <- diagnostics

  class(out) <- "svcLMstack"

  return(out)

}
