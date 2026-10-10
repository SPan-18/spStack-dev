#' Synthetic point-referenced spatial-temporal data
#'
#' @description Dataset of size 500 with spatial coordinates sampled uniformly
#' from the unit square and temporal coordinates sampled uniformly from the
#' unit interval, two covariates, an intercept and a slope of `x1` that vary
#' over space and time, and Gaussian and Poisson responses that share them.
#' @format a \code{data.frame} object with 500 rows and columns
#' \describe{
#'  \item{`s1, s2`}{2-D coordinates in the unit square.}
#'  \item{`t_coords`}{temporal coordinates in the unit interval.}
#'  \item{`x1, x2`}{covariates sampled from the standard normal distribution.}
#'  \item{`y_gauss`}{Gaussian response.}
#'  \item{`y_pois`}{Poisson count response.}
#'  \item{`z1_true`}{true spatial-temporal effect associated with the
#'  intercept.}
#'  \item{`z2_true`}{true spatial-temporal effect associated with `x1`.}
#' }
#' @usage data(simSpaceTime)
#' @details With \eqn{\ell = (s, t)}, the varying intercept is a wave
#' travelling across space over time and the varying slope of \eqn{x_1}
#' changes with \eqn{s_2} and \eqn{t},
#' \deqn{
#' z_1(\ell) = \sin\{2 \pi (s_1 - t)\}, \quad z_2(\ell) = \cos(2 \pi s_2)
#' \cos(\pi t),
#' }
#' deterministic surfaces that are not draws from Gaussian processes. The
#' responses are generated as
#' \deqn{
#' \begin{aligned}
#' y_{\mathrm{gauss}}(\ell) &\sim N(2 + 5 x_1(\ell) - x_2(\ell) + z_1(\ell) +
#' x_1(\ell) z_2(\ell), 0.5^2),\\
#' y_{\mathrm{pois}}(\ell) &\sim \mathrm{Poisson}(\exp\{2 - 0.5 x_1(\ell) +
#' 0.3 x_2(\ell) + z_1(\ell) + x_1(\ell) z_2(\ell)\}).
#' \end{aligned}
#' }
#' This data can be generated with the code given in the example.
#' @seealso [simSpatial]
#' @examples
#' set.seed(1726)
#' n <- 500
#' s1 <- runif(n)
#' s2 <- runif(n)
#' t_coords <- runif(n)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' z1 <- sin(2 * pi * (s1 - t_coords))
#' z2 <- cos(2 * pi * s2) * cos(pi * t_coords)
#' dat <- data.frame(
#'   s1 = s1, s2 = s2, t_coords = t_coords, x1 = x1, x2 = x2,
#'   y_gauss = rnorm(n, 2 + 5 * x1 - x2 + z1 + x1 * z2, sd = 0.5),
#'   y_pois = rpois(n, exp(2 - 0.5 * x1 + 0.3 * x2 + z1 + x1 * z2)),
#'   z1_true = z1, z2_true = z2
#' )
#' all.equal(dat, simSpaceTime)
#' @author Soumyakanti Pan <span18@ucla.edu>
"simSpaceTime"
