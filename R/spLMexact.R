#' Univariate Bayesian spatial linear model
#'
#' @description Fits a Bayesian spatial linear model with spatial process
#'  parameters and the noise-to-spatial variance ratio fixed to a value supplied
#'  by the user. The output contains posterior samples of the fixed effects,
#'  variance parameter, spatial random effects and, if required, leave-one-out
#'  predictive densities.
#' @details Suppose \eqn{\chi = (s_1, \ldots, s_n)} denotes the \eqn{n}
#' spatial locations the response \eqn{y} is observed. With this function, we
#' fit a conjugate Bayesian hierarchical spatial model
#' \deqn{
#' \begin{aligned}
#' y \mid z, \beta, \sigma^2 &\sim N(X\beta + z, \sigma^2 I_n), \quad
#' z \mid \sigma^2_z \sim N(0, \sigma^2_z R(\chi; \phi, \nu)), \\
#' \beta \mid \sigma^2 &\sim N(\mu_\beta, \sigma^2 V_\beta), \quad
#' \sigma^2 \sim \mathrm{IG}(a_\sigma, b_\sigma)
#' \end{aligned}
#' }
#' where we fix the noise-to-spatial variance ratio
#' \eqn{\delta^2 = \sigma^2 / \sigma^2_z}, the spatial process parameters
#' \eqn{\phi} and \eqn{\nu}, and the hyperparameters \eqn{\mu_\beta},
#' \eqn{V_\beta}, \eqn{a_\sigma} and \eqn{b_\sigma}. If \code{priors = "flat"},
#' we instead assign the prior \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}. We utilize
#' a composition sampling strategy to sample the model parameters from their
#' joint posterior distribution which can be written as
#' \deqn{
#' p(\sigma^2, \beta, z \mid y) = p(\sigma^2 \mid y) \times
#' p(\beta \mid \sigma^2, y) \times p(z \mid \beta, \sigma^2, y).
#' }
#' We proceed by first sampling \eqn{\sigma^2} from its marginal posterior,
#' then given the samples of \eqn{\sigma^2}, we sample \eqn{\beta} and
#' subsequently, we sample \eqn{z} conditioned on the posterior samples of
#' \eqn{\beta} and \eqn{\sigma^2} (Banerjee 2020). Posterior samples of the
#' spatial variance are obtained as \eqn{\sigma^2_z = \sigma^2 / \delta^2}.
#' @param formula a symbolic description of the regression model to be fit.
#'  See example below.
#' @param data an optional data frame containing the variables in the model.
#'  If not found in \code{data}, the variables are taken from
#'  \code{environment(formula)}, typically the environment from which
#'  \code{spLMexact} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#'  coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#'  to model the spatial dependence structure among the observations. Supported
#'  covariance model key words are: \code{'exponential'} and \code{'matern'}.
#'  See below for details.
#' @param priors either \code{"flat"} (default), which assigns the prior
#'  \eqn{p(\beta, \sigma^2) \propto 1/\sigma^2}, or a list with tags
#'  \code{beta.norm} (a list containing \eqn{\mu_\beta} and \eqn{V_\beta})
#'  and/or \code{sigma.sq.ig} (a vector containing \eqn{a_\sigma} and
#'  \eqn{b_\sigma}). A component not supplied in the list receives its flat
#'  prior, \eqn{p(\beta) \propto 1} or \eqn{p(\sigma^2) \propto 1/\sigma^2}.
#' @param spParams fixed value of spatial process parameters.
#' @param noise_sp_ratio noise-to-spatial variance ratio.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd logical. If `loopd=TRUE`, returns leave-one-out predictive
#'  densities, using method as given by \code{loopd.method}. Default is
#'  \code{FALSE}.
#' @param loopd.method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#'  valid inputs are `'exact'` and `'PSIS'`. The option `'exact'` corresponds to
#'  exact leave-one-out predictive densities which requires computation almost
#'  equivalent to fitting the model \eqn{n} times. The option `'PSIS'` is
#'  faster and finds approximate leave-one-out predictive densities using
#'  Pareto-smoothed importance sampling (Gelman *et al.* 2024).
#' @param verbose logical. If \code{verbose = TRUE}, prints model description.
#' @param ... currently no additional argument.
#' @return An object of class \code{spLMexact}, which is a list with the
#'  following tags -
#' \describe{
#' \item{samples}{a list of length 4, containing posterior samples of fixed
#'  effects (\code{beta}), measurement error variance (\code{sigmaSq}),
#'  spatial variance (\code{sigmaSq.z}), and spatial effects (\code{z}).}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#'  densities.}
#' \item{model.params}{Values of the fixed parameters that includes
#'  \code{phi} (spatial decay), \code{nu} (spatial smoothness; \code{NA} for
#'  the exponential correlation function) and \code{noise_sp_ratio}
#'  (noise-to-spatial variance ratio).}
#' \item{diagnostics}{a list of fit diagnostics, obtained from quantities the
#'  fit computes anyway. Element \code{numerical} is a data frame with one row
#'  and columns \code{min.pivot} (the smallest relative Cholesky pivot of the
#'  \eqn{n \times n}{n x n} factorizations; values below 1e-8 indicate a
#'  nearly singular covariance matrix), \code{min.cor} and \code{max.cor} (the
#'  correlations of the two farthest-apart and of the two closest locations;
#'  values of \code{min.cor} above 0.95 suggest an effective range far
#'  exceeding the extent of the data, values of \code{max.cor} below 0.05
#'  nearly uncorrelated locations). If \code{loopd.method = 'PSIS'}, element
#'  \code{pareto} is a list with the Pareto \eqn{k} diagnostic values of the
#'  leave-one-out predictive densities (\code{k}), the threshold above which
#'  they are unreliable (\code{threshold}, Vehtari *et al.* 2024) and the
#'  number of values above it (\code{n.high}). If \code{verbose = TRUE}, a
#'  "Diagnostics" section is printed when any threshold is crossed.}
#' }
#' The return object might include additional data used for subsequent
#' prediction and/or model fit evaluation.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [spLMstack()]
#' @references Banerjee S (2020). "Modeling massive spatial datasets using a
#' conjugate Bayesian linear modeling framework." *Spatial Statistics*, **37**,
#' 100417. ISSN 2211-6753. \doi{10.1016/j.spasta.2020.100417}.
#' @references Vehtari A, Simpson D, Gelman A, Yao Y, Gabry J (2024). "Pareto
#'  Smoothed Importance Sampling." *Journal of Machine Learning Research*,
#'  **25**(72), 1-58. URL \url{https://jmlr.org/papers/v25/19-556.html}.
#' @examples
#' data(simSpatial)
#' dat <- simSpatial[1:100, ]
#'
#' # setup prior list
#' muBeta <- c(0, 0, 0)
#' VBeta <- diag(100, 3)
#' sigmaSqIGa <- 2
#' sigmaSqIGb <- 0.1
#' prior_list <- list(beta.norm = list(muBeta, VBeta),
#'                    sigma.sq.ig = c(sigmaSqIGa, sigmaSqIGb))
#'
#' mod1 <- spLMexact(y_gauss ~ x1 + x2, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   priors = prior_list,
#'                   spParams = list(phi = 6, nu = 0.5),
#'                   noise_sp_ratio = 0.5,
#'                   n.samples = 100,
#'                   loopd = TRUE, loopd.method = "exact")
#'
#' post_beta <- mod1$samples$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#'
#' # compare the posterior medians of the spatial effects with the truth
#' cor(apply(mod1$samples$z, 1, median), dat$z_true)
#' @export
spLMexact <- function(formula, data = parent.frame(), coords, cor.fn,
                      priors = "flat",
                      spParams, noise_sp_ratio, n.samples,
                      loopd = FALSE, loopd.method = "exact",
                      verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

  ##### formula #####
  if(missing(formula)){
    stop("error: formula must be specified!")
  }

  if(inherits(formula, "formula")){
    holder <- parseFormula(formula, data)
    y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    X.names <- holder[[3]]
  } else {
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)

  ## storage mode
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ##### coords #####
  if(!is.matrix(coords)){
    stop("error: coords must n-by-2 matrix of xy-coordinate locations")
  }
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or,
    number of rows is different than data used in the model formula")
  }

  check_no_missing(y = y, X = X, coords = coords)
  check_distinct_coords(coords)

  ## distances are computed in C++ from the coordinates
  storage.mode(coords) <- "double"

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("error: cor.fn must be specified")
  }
  if(!cor.fn %in% c("exponential", "matern")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('exponential', 'matern').")
  }

  ##### priors #####
  # priors = "flat" assigns p(beta, sigma.sq) proportional to 1/sigma.sq; if
  # priors is a list, a component not supplied receives its flat prior, i.e.,
  # p(beta) proportional to 1 or p(sigma.sq) proportional to 1/sigma.sq
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.prior <- "flat"
  sigma.sq.IG <- c(0.0, 0.0)

  if(is.character(priors)){

    if(length(priors) != 1 || tolower(priors) != "flat"){
      stop("error: priors must be either 'flat' or a named list with tags
           'beta.norm' and/or 'sigma.sq.ig'.")
    }

  }else if(is.list(priors)){

    if(is.null(names(priors))){
      stop("error: priors must be either 'flat' or a named list with tags
           'beta.norm' and/or 'sigma.sq.ig'.")
    }
    names(priors) <- tolower(names(priors))
    if(any(!names(priors) %in% c("beta.norm", "sigma.sq.ig"))){
      stop("error: invalid tag(s) in priors: '",
           paste(setdiff(names(priors), c("beta.norm", "sigma.sq.ig")),
                 collapse = "', '"),
           "'. Valid tags are 'beta.norm' and 'sigma.sq.ig'.")
    }

    ## Setup prior for beta
    if("beta.norm" %in% names(priors)){
      beta.Norm <- priors[["beta.norm"]]
      if(!is.list(beta.Norm) || length(beta.Norm) != 2){
        stop("error: beta.Norm must be a list of length 2")
      }
      if(length(beta.Norm[[1]]) != p){
        stop(paste("error: beta.Norm[[1]] must be a vector of length, ", p, ".",
                   sep = ""))
      }
      if(length(beta.Norm[[2]]) != p^2){
        stop(paste("error: beta.Norm[[2]] must be a ", p, "x", p,
                   " covariance matrix.", sep = ""))
      }
      check_cov_matrix(beta.Norm[[2]], p, "the prior covariance of beta (beta.norm[[2]])")
      storage.mode(beta.Norm[[1]]) <- "double"
      storage.mode(beta.Norm[[2]]) <- "double"
      beta.prior <- "normal"
    }

    ## Setup prior for sigma.sq
    if("sigma.sq.ig" %in% names(priors)){
      sigma.sq.IG <- priors[["sigma.sq.ig"]]
      if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){
        stop("error: sigma.sq.IG must be a vector of length 2")
      }
      if(any(sigma.sq.IG <= 0)){
        stop("error: sigma.sq.IG must be a positive vector of length 2")
      }
      sigma.sq.prior <- "ig"
    }

  }else{
    stop("error: priors must be either 'flat' or a named list with tags
         'beta.norm' and/or 'sigma.sq.ig'.")
  }

  ## storage mode
  storage.mode(sigma.sq.IG) <- "double"

  ##### spatial process parameters #####
  phi <- 0
  nu <- 0

  if(missing(spParams)){
    stop("spParams (spatial process parameters) must be supplied.")
  }

  names(spParams) <- tolower(names(spParams))

  if(!"phi" %in% names(spParams)){
    stop("phi must be supplied.")
  }
  phi <- spParams[["phi"]]

  if(!is.numeric(phi) || length(phi) != 1){
    stop("phi must be a numeric scalar.")
  }
  if(phi <= 0){
    stop("phi (decay parameter) must be a positive real number.")
  }

  if(cor.fn == "matern"){

    if (!"nu" %in% names(spParams)) {
      stop("nu (smoothness parameter) must be supplied.")
    }
    nu <- spParams[["nu"]]

    if (!is.numeric(nu) || length(nu) != 1) {
      stop("nu must be a numeric scalar.")
    }
    if (nu <= 0) {
      stop("nu (smoothness parameter) must be a positive real number.")
    }

  }

  ## storage mode
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"

  ##### noise-to-spatial variance ratio #####
  deltasq <- 0

  if(missing(noise_sp_ratio)){
    message("noise_sp_ratio not supplied. Using noise_sp_ratio = 1.")
    deltasq = 1
  }else{
    deltasq <- noise_sp_ratio
    if(!is.numeric(deltasq) || length(deltasq) != 1){
      stop("noise_sp_ratio must be a numeric scalar.")
    }
    if(deltasq <= 0){
      stop("noise_sp_ratio must be a positive real number.")
    }
  }

  ## storage mode
  storage.mode(deltasq) <- "double"

  ##### sampling setup #####

  if (missing(n.samples)) {
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ##### Leave-one-out setup #####

  if(loopd){
    loopd.method <- tolower(loopd.method)
    if(!loopd.method %in% c("exact", "psis")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose
           from c('exact', 'PSIS').")
    }
  }else{
    loopd.method <- "none"
  }

  ## sample size check: a flat prior on beta requires n > p for a proper
  ## posterior, and n - 1 > p for the leave-one-out predictive densities
  if(beta.prior == "flat"){
    if(n <= p){
      stop("error: a flat prior on beta requires n > p; supply beta.norm in
           priors.")
    }
    if(loopd && n <= p + 1){
      stop("error: a flat prior on beta requires n - 1 > p for leave-one-out
           predictive densities; supply beta.norm in priors.")
    }
  }

  ##### main function call #####
  ptm <- proc.time()

  if(loopd){
    samps <- .Call(C_spLMexactLOO, y, X, p, n, coords, beta.prior, beta.Norm,
                   sigma.sq.IG, phi, nu, deltasq, cor.fn, n.samples, loopd,
                   loopd.method, verbose)
  }else{
    samps <- .Call(C_spLMexact, y, X, p, n, coords, beta.prior, beta.Norm,
                   sigma.sq.IG, phi, nu, deltasq, cor.fn, n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  if(beta.prior == "normal"){
    beta.Norm.out <- list(mu = beta.Norm[[1]], V = matrix(beta.Norm[[2]], p, p))
  }else{
    beta.Norm.out <- "flat"
  }
  if(sigma.sq.prior == "ig"){
    sigma.sq.IG.out <- sigma.sq.IG
  }else{
    sigma.sq.IG.out <- "flat"
  }
  out$priors <- list(beta.Norm = beta.Norm.out, sigma.sq.IG = sigma.sq.IG.out)
  out$n.samples <- n.samples
  out$samples <- samps[c("beta", "sigmaSq", "sigmaSq.z", "z")]
  if(loopd){
    out$loopd.method <- loopd.method
    out$loopd <- samps[["loopd"]]
  }
  if(cor.fn == 'matern'){
    out$model.params <- list(phi = phi, nu = nu, noise_sp_ratio = deltasq)
  }else{
    out$model.params <- list(phi = phi, nu = NA_real_, noise_sp_ratio = deltasq)
  }
  out$diagnostics <- list(numerical = collect_diagnostics(list(samps)))
  if(loopd && loopd.method == "psis"){
    out$diagnostics$pareto <- pareto_diagnostics(samps[["loopd.pareto_k"]], n.samples)
  }
  out$run.time <- run.time

  if(verbose){
    print_diagnostics(out$diagnostics,
                      pivot.hint = "nearly coincident locations, a very small decay parameter, or a very small noise-to-spatial variance ratio")
  }

  class(out) <- "spLMexact"

  return(out)

}
