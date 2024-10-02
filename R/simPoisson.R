#' Synthetic point-referenced Poisson count data
#'
#' @description Dataset of size 500, with a Poisson distributed response
#' variable indexed by spatial coordinates sampled uniformly from the unit
#' square. The model includes one covariate and spatial random effects induced
#' by a Matérn covariogram.
#' @format a \code{data.frame} object.
#' \describe{
#'  \item{`s1, s2`}{2-D coordinates; latitude and longitude.}
#'  \item{`x1`}{a covariate sampled from the standard normal distribution.}
#'  \item{`y`}{response vector.}
#'  \item{`z_true`}{true spatial random effects that generated the data.}
#' }
#' @usage data(simPoisson)
#' @details With \eqn{n = 500}, the count data is simulated using
#' \deqn{
#' \begin{aligned}
#' y(s_i) &\sim \mathrm{Poisson}(\lambda(s_i)),
#' i = 1, \ldots, n,\\
#' \log \lambda(s_i) &= x(s_i)^\top \beta + z(s_i)
#' \end{aligned}
#' }
#' where the spatial effects \eqn{z \sim N(0, \sigma^2 R)} with \eqn{R} being a
#' \eqn{n \times n} correlation matrix given by the Matérn covariogram
#' \deqn{
#' R(s, s') = \frac{(\phi |s-s'|)^\nu}{\Gamma(\nu) 2^{\nu - 1}}
#' K_\nu(\phi |s-s'|),
#' }
#' where \eqn{\phi} is the spatial decay parameter and \eqn{\nu} the spatial
#' smoothness parameter. We have sampled the data with \eqn{\beta = (2, -0.5)},
#' \eqn{\phi = 5}, \eqn{\nu = 0.5}, and \eqn{\sigma^2 = 0.4}. This data can be
#' generated with the code as given in the example below.
#' @seealso [simGaussian], [simBinom], [simBinary]
#' @examples
#' set.seed(1729)
#' n <- 500
#' beta <- c(2, -0.5)
#' phi0 <- 5
#' nu0 <- 0.5
#' spParams <- c(phi0, nu0)
#' spvar <- 0.4
#' sim1 <- sim_spData(n = n, beta = beta, cor.fn = "matern",
#'                    spParams = spParams, spvar = spvar, deltasq = deltasq,
#'                    family = "poisson")
#'
#' # Plot an interpolated spatial surface of the true random spatial effects
#' plot1 <- surfaceplot(sim1, coords_name = c("s1", "s2"), var_name = "z_true")
#'
#' # Plot the simulated count data
#' library(ggplot2)
#' plot2 <- ggplot(sim1, aes(x = s1, y = s2)) +
#'   geom_point(aes(color = y), alpha = 0.75) +
#'   scale_color_distiller(palette = "RdYlGn", direction = -1,
#'                         label = function(x) sprintf("%.0f", x)) +
#'   guides(alpha = 'none') + theme_bw() +
#'   theme(axis.ticks = element_line(linewidth = 0.25),
#'         panel.background = element_blank(), panel.grid = element_blank(),
#'         legend.title = element_text(size = 10, hjust = 0.25),
#'         legend.box.just = "center", aspect.ratio = 1)
"simPoisson"