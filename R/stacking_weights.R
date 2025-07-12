#' Optimal stacking weights
#'
#' Obtains optimal stacking weights given leave-one-out predictive densities for
#' each candidate model.
#' @param log_loopd an \eqn{n \times M}{n x M} matrix with \eqn{i}{i}-th row
#'  containing the leave-one-out predictive densities for the \eqn{i}{i}-th
#'  data point for the \eqn{M}{M} candidate models.
#' @param solver specifies the solver to use for obtaining optimal weights.
#'  Default is \code{"ECOS"}. Internally calls
#'  [CVXR::psolve()].
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
#'                   parallel = FALSE, solver = "ECOS", verbose = TRUE)
#'
#' loopd_mat <- do.call('cbind', mod1$loopd)
#' w_hat <- get_stacking_weights(loopd_mat)
#' print(round(w_hat$weights, 4))
#' print(w_hat$status)
#' @importFrom CVXR Maximize Problem Variable psolve log_sum_exp Parameter
#' @references Yao Y, Vehtari A, Simpson D, Gelman A (2018). "Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion)." *Bayesian
#' Analysis*, **13**(3), 917-1007. \doi{10.1214/17-BA1091}.
#' @seealso [CVXR::psolve()], [spLMstack()], [spGLMstack()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
get_stacking_weights <- function(log_loopd, solver = "ECOS"){

  # rescale log leave-one-out predictuve densities for numerical stability
  log_loopd_m <- mean(log_loopd)
  log_loopd <- log_loopd - log_loopd_m
  loopd <- exp(log_loopd)
  M <- ncol(loopd)

  # setup CVX optimization problem and constraints
  w <- CVXR::Variable(M)
  obj <- CVXR::Maximize(sum(log(loopd %*% w)))
  constr <- list(sum(w) == 1, w >= 0)
  prob <- CVXR::Problem(objective = obj, constraints = constr)

  # solve the optimization problem using available solvers
  result <- CVXR::psolve(prob, solver = solver)

  # output
  wts <- as.numeric(result$getValue(w))
  return(list(weights = wts, status = result$status))

}