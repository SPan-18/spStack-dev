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
#' @return A list of length 2.
#' \describe{
#'   \item{\code{weights}}{optimal stacking weights as a numeric vector of
#'   length \eqn{M}{M}}
#'   \item{\code{status}}{solver status, returns \code{"optimal"} if solver
#'   succeeded.}
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
#' print(w_hat$status)
#' @importFrom CVXR Maximize Problem Variable psolve Parameter installed_solvers status value
#' @importFrom loo stacking_weights
#' @references Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion)." *Bayesian
#' Analysis*, **13**(3), 917-1007. \doi{10.1214/17-BA1091}.
#' @seealso [CVXR::psolve()], [spLMstack()], [spGLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
get_stacking_weights <- function(log_loopd, solver = NULL, verbose = TRUE){

  # rescale log predictive densities for numerical stability
  log_loopd_m <- mean(log_loopd)
  log_loopd <- log_loopd - log_loopd_m
  loopd <- exp(log_loopd)

  M <- ncol(loopd)

  # CVXR optimization problem
  w <- CVXR::Variable(M)
  expr <- loopd %*% w

  obj <- CVXR::Maximize(sum(log(expr + 1e-12)))
  constr <- list(
    sum(w) == 1,
    w >= 0
  )

  prob <- CVXR::Problem(objective = obj, constraints = constr)

  # solver preference order
  preferred <- c("CLARABEL", "ECOS", "SCS")
  installed <- CVXR::installed_solvers()
  solvers <- intersect(preferred, installed)

  # allow user override
  if (!is.null(solver)) {
    solvers <- unique(c(solver, solvers))
  }

  out <- NULL
  last_error <- NULL

  # try CVXR solvers
  for (s in solvers) {

    result <- tryCatch({
      CVXR::psolve(prob, solver = s, verbose = verbose)
    }, error = function(e){
      last_error <<- conditionMessage(e)
      NULL
    })

    if (!is.null(result)) {

      w_hat <- CVXR::value(w)
      solver_status <- CVXR::status(prob)
      
      w_hat <- as.numeric(w_hat)
      w_hat[!is.finite(w_hat)] <- 0
      w_hat <- pmax(0, w_hat)

      if (sum(w_hat) > 0) {
        w_hat <- w_hat / sum(w_hat)
        out <- list(
          weights = w_hat,
          status = solver_status,
          solver = s
        )
        break
      }
    }
  }

  # fallback to loo if CVXR fails
  if (is.null(out)) {

    message("CVXR solvers failed. Switching to loo::stacking_weights().")

    w_hat <- tryCatch({
      loo::stacking_weights(loopd)
    }, error = function(e){
      stop(
        "Both CVXR and loo stacking failed. Last CVXR error: ",
        last_error
      )
    })

    w_hat <- as.numeric(w_hat)
    w_hat[!is.finite(w_hat)] <- 0
    w_hat <- pmax(0, w_hat)
    w_hat <- w_hat / sum(w_hat)

    out <- list(
      weights = w_hat,
      status = "loo:optimal",
      solver = "loo"
    )
  }

  return(out)
}