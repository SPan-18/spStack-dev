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
#' y \mid z, \beta, \sigma^2 &\sim N(X\beta + z, \delta^2 \sigma^2 I_n), \quad
#' z \mid \sigma^2 \sim N(0, \sigma^2 R(\chi; \phi, \nu)), \\
#' \beta \mid \sigma^2 &\sim N(\mu_\beta, \sigma^2 V_\beta), \quad
#' \sigma^2 \sim \mathrm{IG}(a_\sigma, b_\sigma)
#' \end{aligned}
#' }
#' where we fix the spatial process parameters \eqn{\phi} and \eqn{\nu}, the
#' noise-to-spatial variance ratio \eqn{\delta^2} and the hyperparameters
#' \eqn{\mu_\beta}, \eqn{V_\beta}, \eqn{a_\sigma} and \eqn{b_\sigma}. We utilize
#' a composition sampling strategy to sample the model parameters from their
#' joint posterior distribution which can be written as
#' \deqn{
#' p(\sigma^2, \beta, z \mid y) = p(\sigma^2 \mid y) \times
#' p(\beta \mid \sigma^2, y) \times p(z \mid \beta, \sigma^2, y).
#' }
#' We proceed by first sampling \eqn{\sigma^2} from its marginal posterior,
#' then given the samples of \eqn{\sigma^2}, we sample \eqn{\beta} and
#' subsequently, we sample \eqn{z} conditioned on the posterior samples of
#' \eqn{\beta} and \eqn{\sigma^2} (Banerjee 2020).
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
#' @param priors a list with each tag corresponding to a parameter name and
#'  containing prior details.
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
#' \item{samples}{a list of length 3, containing posterior samples of fixed
#'  effects (\code{beta}), variance parameter (\code{sigmaSq}), spatial effects
#'  (\code{z}).}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#'  densities.}
#' \item{model.params}{Values of the fixed parameters that includes
#'  \code{phi} (spatial decay), \code{nu} (spatial smoothness) and
#'  \code{noise_sp_ratio} (noise-to-spatial variance ratio).}
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
#' # load data
#' data(simGaussian)
#' dat <- simGaussian[1:100, ]
#'
#' # setup prior list
#' muBeta <- c(0, 0)
#' VBeta <- cbind(c(1.0, 0.0), c(0.0, 1.0))
#' sigmaSqIGa <- 2
#' sigmaSqIGb <- 0.1
#' prior_list <- list(beta.norm = list(muBeta, VBeta),
#'                    sigma.sq.ig = c(sigmaSqIGa, sigmaSqIGb))
#'
#' # supply fixed values of model parameters
#' phi0 <- 3
#' nu0 <- 0.75
#' noise.sp.ratio <- 0.8
#'
#' mod1 <- spLMexact(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   priors = prior_list,
#'                   spParams = list(phi = phi0, nu = nu0),
#'                   noise_sp_ratio = noise.sp.ratio,
#'                   n.samples = 100,
#'                   loopd = TRUE, loopd.method = "exact")
#'
#' beta.post <- mod1$samples$beta
#' z.post.median <- apply(mod1$samples$z, 1, median)
#' dat$z.post.median <- z.post.median
#' plot1 <- surfaceplot(dat, coords_name = c("s1", "s2"),
#'                      var_name = "z_true")
#' plot2 <- surfaceplot(dat, coords_name = c("s1", "s2"),
#'                      var_name = "z.post.median")
#' plot1
#' plot2
#' @export
spLMexact <- function(formula, data = parent.frame(), coords, cor.fn, priors,
                      spParams, noise_sp_ratio, n.samples,
                      loopd = FALSE, loopd.method = "exact",
                      verbose = TRUE, ...){

  ##### check for unused args #####
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if (!i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

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

  coords.D <- 0
  coords.D <- iDist(coords)

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("error: cor.fn must be specified")
  }
  if(!cor.fn %in% c("exponential", "matern")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('exponential', 'matern').")
  }

  ##### priors #####
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.IG <- 0

  if(missing(priors)){
    beta.prior <- "normal"
    beta.Norm <- list(rep(0.0, p), diag(100.0, p))
    sigma.sq.IG <- c(2, 2)
  }else{

    names(priors) <- tolower(names(priors))

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
      beta.prior <- "normal"
    }

    ## Setup prior for sigma.sq
    if(!"sigma.sq.ig" %in% names(priors)){
      stop("error: sigma.sq.IG must be specified")
    }
    sigma.sq.IG <- priors[["sigma.sq.ig"]]

    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){
      stop("error: sigma.sq.IG must be a vector of length 2")
    }
    if(any(sigma.sq.IG <= 0)){
      stop("error: sigma.sq.IG must be a positive vector of length 2")
    }

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
    if(missing(loopd.method)){
      stop("loopd.method must be specified")
    }else{
      loopd.method <- tolower(loopd.method)
    }
    if(!loopd.method %in% c("exact", "psis")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose
           from c('exact', 'PSIS').")
    }
  }else{
    loopd.method <- "none"
  }

  ##### main function call #####
  ptm <- proc.time()

  if(loopd){
    samps <- .Call("spLMexactLOO", y, X, p, n, coords.D, beta.prior, beta.Norm,
                   sigma.sq.IG, phi, nu, deltasq, cor.fn, n.samples, loopd,
                   loopd.method, verbose)
  }else{
    samps <- .Call("spLMexact", y, X, p, n, coords.D, beta.prior, beta.Norm,
                   sigma.sq.IG, phi, nu, deltasq, cor.fn, n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$priors <- list(beta.Norm = list(mu = beta.Norm[[1]],
                                      V = matrix(beta.Norm[[2]], p, p)),
                     sigma.sq.IG = sigma.sq.IG)
  out$n.samples <- n.samples
  out$samples <- samps[c("beta", "sigmaSq", "z")]
  if(loopd){
    out$loopd.method <- loopd.method
    out$loopd <- samps[["loopd"]]
  }
  if(cor.fn == 'matern'){
    out$model.params <- list(phi = phi, nu = nu, noise_sp_ratio = deltasq)
  }else{
    out$model.params <- list(phi = phi, noise_sp_ratio = deltasq)
  }
  out$run.time <- run.time

  class(out) <- "spLMexact"

  return(out)

}
