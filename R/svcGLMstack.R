#' Bayesian spatially varying coefficients generalized linear model using
#' predictive stacking
#'
#' @description Fits Bayesian generalized linear models with spatially varying
#' coefficients on a collection of candidate models constructed from candidate
#' values of the spatial process parameters and the boundary adjustment
#' parameter supplied by the user, and combines inference by stacking
#' predictive densities. See Pan, Zhang, Bradley, and Banerjee (2025) for more
#' details.
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially varying coefficients. See
#' examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{svcGLMstack} is called.
#' @param family Specifies the distribution of the response as a member of the
#' exponential family. Supported options are `'poisson'`, `'binomial'` and
#' `'binary'`.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial dependence structure among the observations. Supported
#' covariance model key words are: \code{'exponential'} and \code{'matern'}.
#' @param process.type a quoted keyword specifying the model for the spatial
#' processes: `'independent'`, `'independent.shared'` or `'multivariate'`. See
#' [svcGLMexact()].
#' @param priors (optional) a list with each tag corresponding to a
#' hyperparameter name and containing hyperprior details. Valid tags include
#' `V.beta`, `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. See
#' [svcGLMexact()]. If not supplied, uses defaults.
#' @param candidate.models an object of class `candidateModels` containing a
#' list of candidate models for stacking, each with tags `phi`, `nu` (only if
#' `cor.fn = 'matern'`), vectors of length \eqn{r} if `process.type =
#' 'independent'` (use `list()` entries in [candidateModels()]) and scalars
#' otherwise, and `boundary`, a scalar. See [candidateModels()] for details.
#' @param n.samples number of samples to be drawn from the posterior
#' distribution.
#' @param loopd.controls a list with details on how leave-one-out predictive
#' densities (LOO-PD) are to be calculated, with tags `method` (`'exact'` or
#' `'CV'`), `CV.K`, `nMC` and `CV.update`, as in [stvcGLMstack()]. If sample
#' size is more than 100, then the default is `'CV'` with `CV.K = 10`.
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
#' @return An object of class \code{svcGLMstack}, which is a list including the
#'  following tags -
#' \describe{
#' \item{`samples`}{a list of length equal to total number of candidate models
#'  with each entry corresponding to a list of length 3, containing posterior
#'  samples of fixed effects (\code{beta}), the spatially varying coefficients
#'  (\code{z}), and fine scale variation `xi` for that model.}
#' \item{`loopd`}{a list of length equal to total number of candidate models with
#' each entry containing leave-one-out predictive densities under that
#' particular model.}
#' \item{`n.models`}{number of candidate models that are fit.}
#' \item{`model.params`}{a list with one element per candidate model, each a
#'  named list of its parameters: \code{phi}, \code{nu} (\code{NA} for the
#'  exponential correlation function) and \code{boundary}.}
#' \item{`stacking.summary`}{a matrix with one row per candidate model,
#'  containing its parameters and its optimal stacking weight, for display.}
#' \item{`stacking.weights`}{a numeric vector of length equal to the number of
#' candidate models storing the optimal stacking weights.}
#' \item{`run.time`}{a \code{proc_time} object with runtime details.}
#' \item{`diagnostics`}{a list of diagnostics, as in [stvcGLMstack()]: element
#' \code{numerical} (columns \code{min.pivot}, \code{min.cor} and
#' \code{max.cor}, see [svcGLMexact()]) and element \code{solver}. If
#' \code{verbose = TRUE}, a "Diagnostics" section is printed if there is an
#' issue.}
#' }
#' This object can be further used to recover posterior samples of the scale
#' parameters in the model with [recoverGLMscale()], to make predictions at new
#' locations with [posteriorPredict()], and to sample from the stacked
#' posterior with [stackedSampler()].
#' @details Candidate models that share the process parameters \eqn{(\phi,
#' \nu)} are fitted together: the spatial correlation matrices and all the
#' pre-processing for the leave-one-out predictive densities are computed once
#' and shared across their boundary values, and the groups of such candidates
#' are the units of parallel work.
#' @seealso [svcGLMexact()], [svcLMstack()], [spGLMstack()], [stvcGLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' *Bayesian Analysis*, **In Press**. \doi{10.1214/25-BA1582}.
#' @importFrom rstudioapi isAvailable
#' @importFrom parallel detectCores
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_lapply
#' @examples
#' set.seed(1234)
#' data(simSpatial)
#' dat <- simSpatial[1:100, ]
#'
#' # processes sharing (phi, nu): scalar candidates
#' mod.list <- candidateModels(list(phi = c(4, 8), nu = c(1.5),
#'                                  boundary = c(0.5, 0.75)), "cartesian")
#' mod1 <- svcGLMstack(y_svc_pois ~ x1 + x2 + (x1), data = dat,
#'                     family = "poisson",
#'                     coords = as.matrix(dat[, c("s1", "s2")]),
#'                     cor.fn = "matern", process.type = "independent.shared",
#'                     candidate.models = mod.list,
#'                     loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
#'                     n.samples = 500)
#' post_samps <- stackedSampler(mod1)
#' @export
svcGLMstack <- function(formula, data = parent.frame(), family, coords,
                        cor.fn, process.type, priors, candidate.models,
                        n.samples, loopd.controls, parallel = FALSE,
                        solver = NULL, verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

  ##### family, formula, process type #####
  dd <- vc_glm_data(formula, data, if(missing(family)) stop("Family not specified") else family)
  family <- dd$family
  n <- dd$n
  p <- dd$p
  r <- dd$r
  process.type <- check_vc_glm_process_type(process.type)
  nR <- vc_nR(process.type, r)

  ##### coords #####
  if(missing(coords)){
    stop("coords must be supplied.")
  }
  if(is.data.frame(coords)){
    coords <- as.matrix(coords)
  }
  if(!is.matrix(coords) || ncol(coords) != 2 || nrow(coords) != n){
    stop("coords must be an n x 2 matrix of spatial coordinates, with n = ", n,
         " the number of observations in the model formula.")
  }
  storage.mode(coords) <- "double"
  check_no_missing(y = dd$y, X = dd$X, X_tilde = dd$X_tilde, n.binom = dd$n.binom, coords = coords)
  check_distinct_coords(coords, what = "spatial locations",
                        hint = "Average the observations at a common location, or model them as spatial-temporal data with stvcGLMstack().")

  ##### correlation function #####
  cor.fn <- check_svc_cor_fn(if(missing(cor.fn)) NULL else cor.fn)

  ##### priors #####
  pr <- parse_vc_glm_priors(if(missing(priors)) NULL else priors, p, r, process.type)

  #### candidate models ####
  cands <- svc_candidates(if(missing(candidate.models)) NULL else candidate.models,
                          cor.fn, nR, process.type, "boundary")
  cand_boundary <- vapply(cands, function(x) x[["boundary"]], numeric(1))
  check_vc_boundaries(cand_boundary)

  #### Leave-one-out setup ####
  lc <- parse_vc_glm_loopd_controls(if(missing(loopd.controls)) NULL else loopd.controls, n)

  ##### sampling setup #####
  if(missing(n.samples)){
    stop("n.samples must be specified.")
  }
  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  # candidate models sharing (phi, nu) are fitted in one call, which builds the
  # spatial correlation matrices and all the pre-processing (full data and
  # leave-one-out/CV subsets) once and shares it across their boundary values
  cand_phi <- lapply(cands, function(x) x[["phi"]])
  cand_nu <- lapply(cands, function(x) svc_c_nu(x[["nu"]]))
  cand_groups <- vc_candidate_groups(cand_phi, cand_nu)

  # fits the candidate models in group g (the routine of stvcGLMstack() with a
  # spatial correlation function: no temporal coordinates, and nu in place of
  # phi_t)
  fit_group <- function(g){
    idx <- cand_groups[[g]]
    .Call(C_stvcGLMexactLOOgrid, dd$y, dd$X, dd$X_tilde, n, p, r, family,
          dd$n.binom, coords, NULL, cor.fn, pr$V.beta,
          pr$nu.beta, pr$nu.z, pr$sigmaSq.xi, pr$IW.scale, process.type,
          cand_phi[[idx[1]]], cand_nu[[idx[1]]], cand_boundary[idx],
          n.samples, TRUE, lc$method, lc$CV.K, lc$nMC, lc$CV.update, 0L)
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

  stack_out <- svc_stacking_summary(cands, cor.fn, nR, "boundary", w_hat)
  if(verbose){
    pretty_print_matrix(stack_out, heading = "STACKING WEIGHTS:")
  }

  model_names <- paste("Model", seq_along(cands))
  loopd_list <- lapply(samps, function(x) x[["loopd"]])
  names(loopd_list) <- paste("Model", seq_along(cands), sep = "")

  diagnostics <- list(numerical = collect_diagnostics(samps, model_names))
  diagnostics$solver <- c(list(used = solver_used, status = solver_status), out_solver_details)
  if(verbose){
    print_diagnostics(diagnostics, weights = w_hat,
                      pivot.hint = "nearly coincident locations, a very small decay parameter phi, or a large smoothness nu")
  }

  samps <- lapply(samps, function(x) x[c("beta", "z", "xi")])
  names(samps) <- names(loopd_list)

  out <- list()
  out$y <- dd$y
  out$X <- dd$X
  out$X.names <- dd$X.names
  out$X.svc.names <- dd$X_tilde.names
  out$family <- family
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$n.samples <- n.samples
  # model parameters of each candidate, read by name (e.g. by posteriorPredict())
  out$model.params <- cands
  names(out$model.params) <- model_names
  out$priors <- list(mu.beta = rep(0, p), V.beta = pr$V.beta, nu.beta = pr$nu.beta,
                     nu.z = pr$nu.z, sigmaSq.xi = pr$sigmaSq.xi, IW.scale = pr$IW.scale)
  out$samples <- samps
  out$loopd <- loopd_list
  out$loopd.method <- lc$controls
  out$n.models <- length(cands)
  out$stacking.summary <- stack_out
  out$stacking.weights <- w_hat
  out$run.time <- run.time
  out$diagnostics <- diagnostics

  class(out) <- "svcGLMstack"

  return(out)

}
