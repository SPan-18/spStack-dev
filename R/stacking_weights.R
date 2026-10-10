#' Optimal stacking weights
#'
#' Obtains optimal stacking weights given leave-one-out predictive densities for
#' each candidate model.
#' @param log_loopd an \eqn{n \times M}{n x M} matrix with \eqn{i}{i}-th row
#'  containing the leave-one-out predictive densities for the \eqn{i}{i}-th
#'  data point for the \eqn{M}{M} candidate models.
#' @param solver specifies the solver to use for obtaining optimal weights.
#'  Default is \code{"CLARABEL"}. Internally calls
#'  [CVXR::psolve()].
#' @param verbose if `TRUE`, prints output of optimization routine.
#' @return A list with elements:
#' \describe{
#'   \item{\code{weights}}{optimal stacking weights as a numeric vector of
#'   length \eqn{M}{M} (\code{NA} if no solver succeeded, see Details).}
#'   \item{\code{status}}{solver status, returns \code{"optimal"} if solver
#'   succeeded, and \code{"failed"} if no solver succeeded.}
#'   \item{\code{solver}}{name of the solver used (\code{"none"} if no
#'   solver succeeded).}
#'   \item{\code{details}}{a list with the installed CVXR solvers
#'   (\code{installed}), the requested solver(s) (\code{requested}) and those
#'   of them not installed (\code{missing.requested}), the order in which the
#'   solvers were tried (\code{search.order}), a data frame of the attempts
#'   with the status or error of each solver (\code{attempts}), and whether
#'   the fallback \code{loo::stacking_weights()} was used (\code{fallback}).}
#' }
#' @details The weights maximize the log score of the stacked leave-one-out
#'  predictive densities (Yao *et al.* 2018) over the simplex, using the CVXR
#'  solvers in the order given above. If none of them reaches an optimal
#'  solution, \code{loo::stacking_weights()} is used as a fallback when the
#'  package \pkg{loo} is installed; otherwise the weights are returned as
#'  \code{NA} with status \code{"failed"}.
#' @examples
#' set.seed(1234)
#' data(simGaussian)
#' dat <- simGaussian[1:100, ]
#'
#' cand.mod <- candidateModels(list(phi = c(1.5, 3),
#'                                  nu = c(0.5, 1),
#'                                  noise_sp_ratio = c(1)), "cartesian")
#'
#' mod1 <- spLMstack(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   candidate.models = cand.mod,
#'                   n.samples = 1000, loopd.method = "exact",
#'                   parallel = FALSE, verbose = TRUE)
#'
#' loopd_mat <- do.call('cbind', mod1$loopd)
#' w_hat <- get_stacking_weights(loopd_mat)
#' print(round(w_hat$weights, 4))
#' print(w_hat$solver)
#' print(w_hat$status)
#' @import CVXR
#' @references Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion)." *Bayesian
#' Analysis*, **13**(3), 917-1007. \doi{10.1214/17-BA1091}.
#' @seealso [CVXR::psolve()], [spLMstack()], [spGLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
get_stacking_weights <- function(log_loopd, solver = NULL, verbose = TRUE){

  if (!is.matrix(log_loopd))
    stop("log_loopd must be a matrix")

  if (any(!is.finite(log_loopd)))
    stop("log_loopd contains non-finite values")

  OPTIMAL_STATUSES <- c("optimal", "optimal_inaccurate")

  ## -----------------------------
  ## Numerical stabilization
  ## -----------------------------

  shift <- max(log_loopd, na.rm = TRUE)
  loopd <- exp(log_loopd - shift)

  M <- ncol(loopd)

  ## -----------------------------
  ## Optimization problem
  ## -----------------------------

  w <- CVXR::Variable(M)
  expr <- loopd %*% w

  obj <- CVXR::Maximize(sum(log(expr + 1e-12)))

  constr <- list(
    sum(w) == 1,
    w >= 0
  )

  prob <- CVXR::Problem(objective = obj, constraints = constr)

  ## -----------------------------
  ## Solver selection and diagnostics
  ## -----------------------------

  installed <- CVXR::installed_solvers()

  if (length(installed) == 0) {
    stop("No CVXR solvers are installed.")
  }

  ## Determine solver order
  missing <- character(0)
  requested_input <- if (is.null(solver)) "DEFAULT (CLARABEL -> ECOS -> SCS)" else solver
  if (!is.null(solver)) {

    missing <- setdiff(solver, installed)

    if (length(missing) > 0) {
      message(
        "Requested solver(s) not installed: ",
        paste(missing, collapse = ", "),
        "\nFalling back to default solver order: CLARABEL -> ECOS -> SCS."
      )
      solver <- NULL
    }
  }

  if (!is.null(solver)) {
    requested <- solver
    solvers <- solver
  } else {
    requested <- "DEFAULT (CLARABEL -> ECOS -> SCS)"
    preferred <- c("CLARABEL", "ECOS", "SCS")
    solvers <- intersect(preferred, installed)

    if (length(solvers) == 0) {
      if (verbose)
        message("No preferred solvers installed; using all available CVXR solvers.")
      solvers <- installed
    }
  }

  ## Diagnostics
  if (verbose) {
    message("--------------------------------------------------")
    message("Solver diagnostics:")
    message("Installed solvers: ", paste(installed, collapse = ", "))
    message("Requested solver: ", paste(requested, collapse = ", "))
    message("Solver search order: ", paste(solvers, collapse = " -> "))
    message("--------------------------------------------------")
  }

  ## -----------------------------
  ## Solver cascade
  ## -----------------------------

  out <- NULL
  last_error <- NULL
  attempts <- data.frame(solver = character(0), status = character(0),
                         stringsAsFactors = FALSE)

  for (s in solvers) {

    result <- tryCatch(
      CVXR::psolve(prob, solver = s, verbose = verbose),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(result)){
      attempts[nrow(attempts) + 1, ] <- c(s, paste0("error: ", last_error))
      next
    }

    ## Solution extraction
    w_hat <- tryCatch(
      as.numeric(CVXR::value(w)),
      error = function(e) NULL
    )

    w_hat[!is.finite(w_hat)] <- 0
    w_hat <- pmax(0, w_hat)

    w_hat_sum <- sum(w_hat)
    if (!is.finite(w_hat_sum) || w_hat_sum <= 0){
      attempts[nrow(attempts) + 1, ] <- c(s, "no valid solution")
      next
    }

    w_hat <- w_hat / w_hat_sum

    ## Status extraction
    solver_status <- .get_cvxr_status(prob, result)

    if (is.na(solver_status))
      solver_status <- "unknown"
    attempts[nrow(attempts) + 1, ] <- c(s, solver_status)

    if (solver_status %in% OPTIMAL_STATUSES) {

      out <- list(
        weights = w_hat,
        status = solver_status,
        solver = paste0("CVXR:", s)
      )

      break
    }
  }

  ## -----------------------------
  ## Fallback solver
  ## -----------------------------

  # loo (in Suggests) is used only here; without it the weights are NA
  if (is.null(out) && !requireNamespace("loo", quietly = TRUE)) {

    if(verbose){
      message("CVXR solvers failed or did not reach optimality, and the",
              " fallback loo::stacking_weights() needs the 'loo' package.")
    }

    out <- list(
      weights = rep(NA_real_, M),
      status = "failed",
      solver = "none"
    )
  }

  if (is.null(out)) {

    if(verbose){
      message("CVXR solvers failed or did not reach optimality.")
      message("Switching to loo::stacking_weights().")
    }

    # loo::stacking_weights() takes the pointwise log predictive densities
    w_hat <- tryCatch(
      loo::stacking_weights(log_loopd),
      error = function(e) {
        stop(
          "Both CVXR and loo stacking failed. Last CVXR error: ",
          last_error, "; loo error: ", conditionMessage(e)
        )
      }
    )

    w_hat <- as.numeric(w_hat)
    w_hat[!is.finite(w_hat)] <- 0

    if (sum(w_hat) > 0)
      w_hat <- w_hat / sum(w_hat)

    out <- list(
      weights = w_hat,
      status = "optimal",
      solver = "loo"
    )
  }

  out$details <- list(installed = installed, requested = requested_input,
                      missing.requested = missing, search.order = solvers,
                      attempts = attempts, fallback = identical(out$solver, "loo"))

  out
}

# ------------------------------------------------------------------
# Internal helper: extract solver status safely across CVXR versions
# ------------------------------------------------------------------
#' @importFrom utils getFromNamespace
.get_cvxr_status <- function(prob, result) {

  status_fun <- try(getFromNamespace("status", "CVXR"), silent = TRUE)

  if (!inherits(status_fun, "try-error")) {
    status_val <- try(status_fun(prob), silent = TRUE)

    if (!inherits(status_val, "try-error") && !is.null(status_val)) {
      return(tolower(as.character(status_val)))
    }
  }

  if (is.list(result) && !is.null(result[["status"]])) {
    return(tolower(as.character(result[["status"]])))
  }

  NA_character_
}