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
#'   length \eqn{M}{M}}
#'   \item{\code{status}}{solver status, returns \code{"optimal"} if solver
#'   succeeded.}
#'   \item{\code{solver}}{name of the solver used.}
#' }
#' @examples
#' set.seed(1234)
#' data(simGaussian)
#' dat <- simGaussian[1:100, ]
#'
#' mod1 <- spLMstack(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   params.list = list(phi = c(1.5, 3),
#'                                      nu = c(0.5, 1),
#'                                      noise_sp_ratio = c(1)),
#'                   n.samples = 1000, loopd.method = "exact",
#'                   parallel = FALSE, verbose = TRUE)
#'
#' loopd_mat <- do.call('cbind', mod1$loopd)
#' w_hat <- get_stacking_weights(loopd_mat)
#' print(round(w_hat$weights, 4))
#' print(w_hat$solver)
#' print(w_hat$status)
#' @import CVXR
#' @importFrom loo stacking_weights
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

  for (s in solvers) {

    result <- tryCatch(
      CVXR::psolve(prob, solver = s, verbose = verbose),
      error = function(e) {
        last_error <- conditionMessage(e)
        NULL
      }
    )

    if (is.null(result))
      next

    ## Solution extraction
    w_hat <- tryCatch(
      as.numeric(CVXR::value(w)),
      error = function(e) NULL
    )

    w_hat[!is.finite(w_hat)] <- 0
    w_hat <- pmax(0, w_hat)

    w_hat_sum <- sum(w_hat)
    if (!is.finite(w_hat_sum) || w_hat_sum <= 0)
      next

    w_hat <- w_hat / w_hat_sum

    ## Status extraction
    solver_status <- .get_cvxr_status(prob, result)

    if (is.na(solver_status))
      solver_status <- "unknown"

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

  if (is.null(out)) {

    if(verbose){
      message("CVXR solvers failed or did not reach optimality.")
      message("Switching to loo::stacking_weights().")
    }

    w_hat <- tryCatch(
      loo::stacking_weights(loopd),
      error = function(e) {
        stop(
          "Both CVXR and loo stacking failed. Last CVXR error: ",
          last_error
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