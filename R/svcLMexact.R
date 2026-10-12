#' Bayesian spatially varying coefficients linear model
#'
#' @description Fits a Bayesian linear model with spatially varying
#' coefficients for a Gaussian response, with the spatial process parameters
#' and the noise-to-spatial variance ratios fixed to values supplied by the
#' user. The output contains exact posterior samples of the fixed effects, the
#' noise and process variances, the spatially varying coefficients and, if
#' required, leave-one-out predictive densities.
#' @details Suppose \eqn{\chi = (s_1, \ldots, s_n)} denotes the \eqn{n}
#' spatial locations in \eqn{\mathcal{S} \subset \mathbb{R}^2} at which the
#' response \eqn{y} is observed. With this function, we fit the conjugate
#' Bayesian hierarchical model
#' \deqn{
#' \begin{aligned}
#' y(s) &= x(s)^\top \beta + \tilde{x}(s)^\top z(s) +
#' \epsilon(s), \quad \epsilon(s) \sim N(0, \sigma^2),\\
#' z_j &\mid \sigma^2 \sim N(0, \sigma^2_{z_j} R(\chi; \phi_j, \nu_j)),
#' \quad \sigma^2_{z_j} = \sigma^2 / \delta^2_j, \quad j = 1, \ldots, r,\\
#' \beta &\mid \sigma^2 \sim N(\mu_\beta, \sigma^2 V_\beta), \quad
#' \sigma^2 \sim \mathrm{IG}(a_\sigma, b_\sigma),
#' \end{aligned}
#' }
#' where \eqn{\tilde{x}(s)} denotes the \eqn{r} covariates with spatially
#' varying coefficients, the processes \eqn{z_1, \ldots, z_r} are independent,
#' and \eqn{R(\chi; \phi, \nu)} is the exponential or Matern spatial
#' correlation matrix with decay \eqn{\phi} and smoothness \eqn{\nu}. We fix the
#' noise-to-spatial variance ratios \eqn{\delta^2_j = \sigma^2 /
#' \sigma^2_{z_j}}, the process parameters \eqn{\phi_j} and \eqn{\nu_j}, and the
#' hyperparameters \eqn{\mu_\beta}, \eqn{V_\beta}, \eqn{a_\sigma} and
#' \eqn{b_\sigma}. If \code{process.type = 'independent'}, each process has its
#' own \eqn{(\phi_j, \nu_j, \delta^2_j)}; if \code{process.type =
#' 'independent.shared'}, they share one. If \code{priors = "flat"}, we instead
#' assign the prior \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}.
#'
#' The joint posterior distribution is available in closed form and is sampled
#' exactly by composition,
#' \deqn{
#' p(\sigma^2, \beta, z \mid y) = p(\sigma^2 \mid y) \times
#' p(\beta \mid \sigma^2, y) \times p(z \mid \beta, \sigma^2, y),
#' }
#' where \eqn{\sigma^2 \mid y} is inverse-gamma and the other two are Gaussian.
#' All of them depend on the data through the \eqn{n \times n}{n x n} matrix
#' \eqn{V_y = I_n + \sum_j \delta_j^{-2} D_j R_j D_j}, with \eqn{D_j =
#' \mathrm{diag}(\tilde{x}_j)}. The \eqn{nr}-dimensional vector \eqn{z} is
#' drawn by a prior draw followed by a kriging correction (Matheron's rule;
#' Bhattacharya, Chakraborty and Mallick 2016), which needs only \eqn{n \times
#' n}{n x n} Cholesky factorizations. The exact leave-one-out predictive
#' densities are obtained in closed form from the same factorizations. The
#' computations are those of [stvcLMexact()], with a spatial correlation
#' function in place of the spatial-temporal one.
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially varying coefficients. See
#' examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{svcLMexact} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial dependence structure among the observations. Supported
#' covariance model key words are: \code{'exponential'} and \code{'matern'}.
#' See below for details.
#' @param process.type a quoted keyword specifying the model for the spatial
#' processes of the varying coefficients. Supported keywords are
#' `'independent'`, independent processes with their own process parameters
#' and noise-to-spatial variance ratios, and `'independent.shared'`,
#' independent processes that share common process parameters and a common
#' noise-to-spatial variance ratio.
#' @param spParams fixed values of the spatial process parameters, a list with
#' tags `phi` (spatial decay) and, if `cor.fn = 'matern'`, `nu` (smoothness).
#' If `process.type = 'independent'`, each is a vector of length \eqn{r},
#' otherwise a scalar.
#' @param noise_sp_ratio noise-to-spatial variance ratio(s) \eqn{\delta^2_j}: a
#' vector of length \eqn{r} if `process.type = 'independent'`, otherwise a
#' scalar. Default is 1.
#' @param priors either \code{"flat"} (default), which assigns the prior
#' \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}, or a list with tags
#' \code{beta.norm} (a list containing \eqn{\mu_\beta} and \eqn{V_\beta})
#' and/or \code{sigma.sq.ig} (a vector containing \eqn{a_\sigma} and
#' \eqn{b_\sigma}). A component not supplied in the list receives its flat
#' prior, \eqn{p(\beta) \propto 1} or \eqn{p(\sigma^2) \propto 1/\sigma^2}.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd logical. If `loopd=TRUE`, returns leave-one-out predictive
#' densities, using method as given by \code{loopd.method}. Default is
#' \code{FALSE}.
#' @param loopd.method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#' valid inputs are `'exact'` and `'PSIS'`. The option `'exact'` finds the exact
#' leave-one-out predictive densities in closed form, at the cost of about one
#' additional \eqn{n \times n}{n x n} triangular inversion. The option `'PSIS'`
#' finds approximate leave-one-out predictive densities using Pareto-smoothed
#' importance sampling (Vehtari *et al.* 2024); with many latent effects
#' (\eqn{nr}), its Pareto \eqn{k} diagnostics are often high and `'exact'` is
#' recommended.
#' @param verbose logical. If \code{verbose = TRUE}, prints model description.
#' @param ... currently no additional argument.
#' @return An object of class \code{svcLMexact}, which is a list with the
#' following tags -
#' \describe{
#' \item{samples}{a list of length 4, containing posterior samples of fixed
#' effects (\code{beta}, a \eqn{p \times} \code{n.samples} matrix), the noise
#' variance (\code{sigmaSq}), the process variances (\code{sigmaSq.z}, an
#' \eqn{r \times} \code{n.samples} matrix if \code{process.type =
#' 'independent'}, otherwise a vector), and the spatially varying coefficients
#' (\code{z}, an \eqn{nr \times} \code{n.samples} matrix whose rows
#' \eqn{(j-1)n + 1, \ldots, jn} correspond to the \eqn{j}-th varying
#' coefficient).}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#' densities.}
#' \item{model.params}{Values of the fixed parameters: \code{phi}, \code{nu}
#' (\code{NA} for the exponential correlation function) and
#' \code{noise_sp_ratio}.}
#' \item{diagnostics}{a list of fit diagnostics, obtained from quantities the
#' fit computes anyway. Element \code{numerical} is a data frame with one row
#' (one row per process, if \code{process.type = 'independent'}) and columns
#' \code{min.pivot} (the smallest relative Cholesky pivot of the spatial
#' correlation matrix and of \eqn{V_y}; values below 1e-8 indicate a nearly
#' singular matrix), \code{min.cor} and \code{max.cor} (the correlations of the
#' two farthest-apart and of the two closest locations; values of
#' \code{min.cor} above 0.95 suggest an effective range far exceeding the
#' extent of the data, values of \code{max.cor} below 0.05 nearly uncorrelated
#' locations). If \code{loopd.method = 'PSIS'}, element \code{pareto} is a list
#' with the Pareto \eqn{k} diagnostic values (\code{k}), the threshold above
#' which they are unreliable (\code{threshold}) and the number of values above
#' it (\code{n.high}). If \code{verbose = TRUE}, a "Diagnostics" section is
#' printed when any threshold is crossed.}
#' }
#' The return object might include additional data used for subsequent
#' prediction and/or model fit evaluation.
#' @seealso [svcLMstack()], [svcGLMexact()], [spLMexact()], [stvcLMexact()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @references Bhattacharya A, Chakraborty A, Mallick BK (2016). "Fast sampling
#' with Gaussian scale mixture priors in high-dimensional regression."
#' *Biometrika*, **103**(4), 985-991. \doi{10.1093/biomet/asw042}.
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto
#'  Smoothed Importance Sampling." *Journal of Machine Learning Research*,
#'  **25**(72), 1-58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @examples
#' data(simSpatial)
#' dat <- simSpatial[1:200, ]
#'
#' # varying intercept and varying slope of x1, each with its own parameters
#' mod1 <- svcLMexact(y_svc ~ x1 + x2 + (x1), data = dat,
#'                    coords = as.matrix(dat[, c("s1", "s2")]),
#'                    cor.fn = "matern", process.type = "independent",
#'                    spParams = list(phi = c(8, 4), nu = c(1.5, 1.5)),
#'                    noise_sp_ratio = c(0.25, 0.25),
#'                    n.samples = 500, loopd = TRUE, verbose = FALSE)
#'
#' # rows 1, ..., n of z hold the varying intercept, rows n+1, ..., 2n the
#' # varying slope of x1
#' n <- nrow(dat)
#' z_hat <- apply(mod1$samples$z, 1, median)
#' cor(z_hat[1:n], dat$z2_true)
#' cor(z_hat[n + 1:n], dat$z_true)
#' @export
svcLMexact <- function(formula, data = parent.frame(), coords, cor.fn,
                       process.type, spParams, noise_sp_ratio,
                       priors = "flat", n.samples, loopd = FALSE,
                       loopd.method = "exact", verbose = TRUE, ...){

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

  ##### correlation function and spatial process parameters #####
  cor.fn <- check_svc_cor_fn(if(missing(cor.fn)) NULL else cor.fn)
  sp <- svc_sp_params(if(missing(spParams)) NULL else spParams, cor.fn, nR, process.type)

  ##### priors #####
  pr <- parse_lm_priors(priors, p)

  ##### noise-to-spatial variance ratios #####
  if(missing(noise_sp_ratio)){
    message("noise_sp_ratio not supplied. Using noise_sp_ratio = 1.")
    noise_sp_ratio <- rep(1, nR)
  }
  deltasq <- check_vc_param(noise_sp_ratio, "noise_sp_ratio", nR, process.type)

  ##### sampling setup #####
  if(missing(n.samples)){
    stop("n.samples must be specified.")
  }
  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ##### Leave-one-out setup #####
  if(loopd){
    loopd.method <- tolower(loopd.method)
    if(!loopd.method %in% c("exact", "psis")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose from c('exact', 'PSIS').")
    }
  }else{
    loopd.method <- "none"
  }
  storage.mode(loopd) <- "integer"

  ## sample size check: a flat prior on beta requires n > p for a proper
  ## posterior, and n - 1 > p for the leave-one-out predictive densities
  if(pr$beta.prior == "flat"){
    if(n <= p){
      stop("a flat prior on beta requires n > p; supply beta.norm in priors.")
    }
    if(loopd && n <= p + 1){
      stop("a flat prior on beta requires n - 1 > p for leave-one-out predictive densities; supply beta.norm in priors.")
    }
  }

  ##### main function call #####
  ptm <- proc.time()

  # the routine of stvcLMexact() with a spatial correlation function: no
  # temporal coordinates, and nu in place of phi_t
  samps <- .Call(C_stvcLMexact, dd$y, dd$X, dd$X_tilde, n, p, r,
                 dd$sp_coords, NULL, cor.fn, process.type,
                 sp$phi, svc_c_nu(sp$nu), pr$beta.prior, pr$beta.Norm, pr$sigma.sq.IG,
                 deltasq, n.samples, loopd, loopd.method, verbose)

  run.time <- proc.time() - ptm

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
  out$samples <- samps[c("beta", "sigmaSq", "sigmaSq.z", "z")]
  if(loopd){
    out$loopd.method <- loopd.method
    out$loopd <- samps[["loopd"]]
  }
  out$model.params <- list(phi = sp$phi, nu = sp$nu, noise_sp_ratio = deltasq)
  out$diagnostics <- list(numerical = collect_diagnostics(list(samps)))
  if(loopd && loopd.method == "psis"){
    out$diagnostics$pareto <- pareto_diagnostics(samps[["loopd.pareto_k"]], n.samples)
  }
  out$run.time <- run.time

  if(verbose){
    print_diagnostics(out$diagnostics,
                      pivot.hint = "nearly coincident locations, very small decay parameters phi, or very small noise-to-spatial variance ratios")
  }

  class(out) <- "svcLMexact"

  return(out)

}
