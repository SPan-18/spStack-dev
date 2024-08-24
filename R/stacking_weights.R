#' Optimal stacking weights
#'
#' Obtains optimal stacking weights given leave-one-out predictive densities for
#' each candidate model.
#' @param log_loopd an \eqn{n \times M}{n x M} matrix with \eqn{i}{i}-th row
#'  containing the leave-one-out predictive densities for the \eqn{i}{i}-th
#'  data point for the \eqn{M}{M} candidate models.
#' @param solver specifies the solver to use for obtaining optimal weights.
#'  Default is \code{"ECOS"}. Internally calls
#'  [CVXR::solve()].
#' @return A list of length 2.
#' \describe{
#'   \item{\code{weights}}{optimal stacking weights as a numeric vecor of length
#'   \eqn{M}{M}}
#'   \item{\code{status}}{solver status. \code{optimal} if solver succeeded,
#'   otherwise throws an error.}
#' }
#' @importFrom CVXR Maximize Problem Variable psolve log_sum_exp Parameter
#' @references Yao Y, Vehtari A, Simpson D, Gelman A (2018). “Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion).” *Bayesian
#' Analysis*, **13**(3), 917 – 1007. \url{https://doi.org/10.1214/17-BA1091}.
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

# stacking_weights <- function(log_loopd, solver){

#   # rescale log leave-one-out predictive densities for numerical stability
#   log_loopd_m <- mean(log_loopd)
#   log_loopd <- log_loopd - log_loopd_m
#   # loopd <- exp(log_loopd)
#   M <- ncol(log_loopd)
#   n <- nrow(log_loopd)

#   # setup CVX optimization problem and constraints
#   w <- CVXR::Variable(M)
#   log_w <- log(w + 1e-8)

#   log_sum_exp_list <- list()

#   for(i in 1:n){
#     log_sum_exp_list[[i]] <- log_sum_exp(log_loopd[i, ] + log_w)
#   }

#   log_sum_exp <- do.call("rbind", log_sum_exp_list)

#   # Sum the resulting values into a single CVXR expression
#   objective_expr <- sum(log_sum_exp)

#   obj <- CVXR::Minimize(objective_expr)

#   constr <- list(w >= 0, sum(w) == 1)
#   prob <- CVXR::Problem(objective = obj, constraints = constr)

#   # solve the optimization problem using available solvers
#   result <- CVXR::psolve(prob, solver = solver)

#   # output
#   wts <- as.numeric(result$getValue(w))
#   return(list(weights = wts, status = result$status))

# }