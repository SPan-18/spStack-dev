#' Bayesian spatially-temporally varying generalized linear model
#'
#' @description Fits a Bayesian generalized linear model with
#' spatially-temporally varying coefficients under fixed values of
#' spatial-temporal process parameters and some auxiliary model parameters. The
#' output contains posterior samples of the fixed effects, spatial-temporal
#' random effects and, if required, finds leave-one-out predictive densities.
#' @details With this function, we fit a Bayesian hierarchical
#' spatially-temporally varying generalized linear model by sampling exactly
#' from the joint posterior distribution utilizing the generalized conjugate
#' multivariate distribution theory (Bradley and Clinch 2024). Suppose
#' \eqn{\chi = (\ell_1, \ldots, \ell_n)} denotes the \eqn{n} spatial-temporal
#' co-ordinates in \eqn{\mathcal{L} = \mathcal{S} \times \mathcal{T}}, the
#' response \eqn{y} is observed. Let \eqn{y(\ell)} be the outcome at the
#' co-ordinate \eqn{\ell} endowed with a probability law from the natural
#' exponential family, which we denote by
#' \deqn{
#' y(\ell) \sim \mathrm{EF}(x(\ell)^\top \beta + \tilde{x}(\ell)^\top z(\ell);
#' b(\ell), \psi)
#' }
#' for some positive parameter \eqn{b(\ell) > 0} and unit log partition function
#' \eqn{\psi}. Here, \eqn{\tilde{x}(\ell)} denotes covariates with
#' spatially-temporally varying coefficients We consider the following response
#' models based on the input supplied to the argument `family`.
#' \describe{
#' \item{`'poisson'`}{It considers point-referenced Poisson responses
#' \eqn{y(\ell) \sim \mathrm{Poisson}(e^{x(\ell)^\top \beta +
#' \tilde{x}(\ell)^\top z(\ell)})}. Here, \eqn{b(\ell) = 1} and
#' \eqn{\psi(t) = e^t}.}
#' \item{`'binomial'`}{It considers point-referenced binomial counts
#' \eqn{y(\ell) \sim \mathrm{Binomial}(m(\ell), \pi(\ell))} where, \eqn{m(\ell)}
#' denotes the total number of trials and probability of success
#' \eqn{\pi(\ell) = \mathrm{ilogit}(x(\ell)^\top \beta + \tilde{x}(\ell)^\top
#' z(\ell))} at spatial-temporal co-ordinate \eqn{\ell}. Here, \eqn{b = m(\ell)}
#' and \eqn{\psi(t) = \log(1+e^t)}.}
#' \item{`'binary'`}{It considers point-referenced binary data (0 or, 1) i.e.,
#' \eqn{y(\ell) \sim \mathrm{Bernoulli}(\pi(\ell))}, where probability of
#' success \eqn{\pi(\ell) = \mathrm{ilogit}(x(\ell)^\top \beta +
#' \tilde{x}(\ell)^\top z(\ell))} at spatial-temporal co-ordinate \eqn{\ell}.
#' Here, \eqn{b(\ell) = 1} and \eqn{\psi(t) = \log(1 + e^t)}.}
#' }
#' The hierarchical model is given as
#' \deqn{
#' \begin{aligned}
#' y(\ell_i) &\mid \beta, z, \xi \sim EF(x(\ell_i)^\top \beta +
#' \tilde{x}(\ell_i)^\top z(s_i) + \xi_i - \mu_i; b_i, \psi_y),
#' i = 1, \ldots, n\\
#' \xi &\mid \beta, z, \sigma^2_\xi, \alpha_\epsilon \sim
#' \mathrm{GCM_c}(\cdots),\\
#' \beta &\mid \sigma^2_\beta \sim N(0, \sigma^2_\beta V_\beta), \quad
#' \sigma^2_\beta \sim \mathrm{IG}(\nu_\beta/2, \nu_\beta/2)\\
#' z_j &\mid \sigma^2_{z_j} \sim N(0, \sigma^2_{z_j} R(\chi; \phi_s, \phi_t)),
#' \quad \sigma^2_{z_j} \sim \mathrm{IG}(\nu_z/2, \nu_z/2), j = 1, \ldots, r
#' \end{aligned}
#' }
#' where \eqn{\mu = (\mu_1, \ldots, \mu_n)^\top} denotes the discrepancy
#' parameter. We fix the spatial-temporal process parameters \eqn{\phi_s} and
#' \eqn{\phi_t} and the hyperparameters \eqn{V_\beta}, \eqn{\nu_\beta},
#' \eqn{\nu_z} and \eqn{\sigma^2_\xi}. The term \eqn{\xi} is known as the
#' fine-scale variation term which is given a conditional generalized conjugate
#' multivariate distribution as prior. For more details, see Pan *et al.* 2024.
#' Default values for \eqn{V_\beta}, \eqn{\nu_\beta}, \eqn{\nu_z},
#' \eqn{\sigma^2_\xi} are diagonal with each diagonal element 100, 2.1, 2.1 and
#' 0.1 respectively.
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially-temporally varying
#' coefficients. See examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{stvcGLMexact} is called.
#' @param family Specifies the distribution of the response as a member of the
#' exponential family. Supported options are `'poisson'`, `'binomial'` and
#' `'binary'`.
#' @param sp_coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' spatial coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param time_coords an \eqn{n \times 1}{n x 1} matrix of the observation
#' temporal coordinates in \eqn{\mathcal{T} \subseteq [0, \infty)}.
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial-temporal dependence structure among the observations.
#' Supported covariance model key words are: \code{'gneiting-decay'} (Gneiting
#' and Guttorp 2010). See below for details.
#' @param priors (optional) a list with each tag corresponding to a
#' hyperparameter name and containing hyperprior details. Valid tags include
#' `V.beta`, `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
#' and `nu.z` must be at least 2.1. If not supplied, uses defaults.
#' @param process.type a quoted keyword specifying the model for the
#' spatial-temporal process. Supported keywords are `'independent'` which
#' indicates independent processes for each varying coefficients characterized
#' by different process parameters, `independent.shared` implies independent
#' processes for the varying coefficients that shares common process parameters,
#' and `multivariate` implies correlated processes for the varying coefficients
#' modeled by a multivariate Gaussian process with an inverse-Wishart prior on
#' the correlation matrix. The input for `sptParams` and `priors` must be given
#' accordingly.
#' @param sptParams fixed values of spatial-temporal process parameters in
#' usually a list of length 2. If `cor.fn='gneiting-decay'`, then it is a list
#' of length 2 with tags `phi_s` and `phi_t`. If `process.type='independent'`,
#' then `phi_s` and `phi_t` contain fixed values of the \eqn{r} spatial-temporal
#' processes, otherwise they will contain scalars. See examples below.
#' @param boundary Specifies the boundary adjustment parameter. Must be a real
#' number between 0 and 1. Default is 0.5.
#' @param n.samples number of posterior samples to be generated.
#' @param loopd logical. If `loopd=TRUE`, returns leave-one-out predictive
#' densities, using method as given by \code{loopd.method}. Default is
#' \code{FALSE}.
#' @param loopd.method character. Ignored if `loopd=FALSE`. If `loopd=TRUE`,
#' valid inputs are `'exact'`, `'CV'`. The option `'exact'` corresponds to exact
#' leave-one-out predictive densities which requires computation almost
#' equivalent to fitting the model \eqn{n} times. The options `'CV'` is faster
#' as it implements \eqn{K}-fold cross validation to find approximate
#' leave-one-out predictive densities (Vehtari *et al.* 2017).
#' @param CV.K An integer between 10 and 20. Considered only if
#' `loopd.method='CV'`. Default is 10 (as recommended in Vehtari *et. al* 2017).
#' @param loopd.nMC Number of Monte Carlo samples to be used to evaluate
#' leave-one-out predictive densities when `loopd.method` is set to either
#' 'exact' or 'CV'.
#' @param verbose logical. If \code{verbose = TRUE}, prints model description.
#' @param ... currently no additional argument.
#' @return An object of class \code{stvcGLMexact}, which is a list with the
#'  following tags -
#' \describe{
#' \item{priors}{details of the priors used, containing the values of the
#' boundary adjustment parameter (`boundary`), the variance parameter of the
#' fine-scale variation term (`simasq.xi`) and others.}
#' \item{samples}{a list of length 3, containing posterior samples of fixed
#'  effects (\code{beta}), spatial-temporal effects (\code{z}) and the
#' fine-scale variation term (\code{xi}). The element with tag \code{z} will
#' again be a list of length \eqn{r}, each containing posterior samples of the
#' spatial-temporal random effects corresponding to each varying coefficient.}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#' densities.}
#' \item{model.params}{Values of the fixed parameters that includes \code{phi_s}
#' (spatial decay), \code{phi_t} (temporal smoothness).}
#' }
#' The return object might include additional data that can be used for
#' subsequent prediction and/or model fit evaluation.
#' @seealso [spGLMexact()]
#' @author Soumyakanti Pan <span18@ucla.edu>
#' @references Bradley JR, Clinch M (2024). "Generating Independent Replicates
#' Directly from the Posterior Distribution for a Class of Spatial Hierarchical
#' Models." *Journal of Computational and Graphical Statistics*, **0**(0), 1-17.
#' \doi{10.1080/10618600.2024.2365728}.
#' @references T. Gneiting and P. Guttorp (2010). "Continuous-parameter
#' spatio-temporal processes." In *A.E. Gelfand, P.J. Diggle, M. Fuentes, and
#' P Guttorp, editors, Handbook of Spatial Statistics*, Chapman & Hall CRC
#' Handbooks of Modern Statistical Methods, p 427–436. Taylor and Francis.
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2024). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' \doi{10.48550/arXiv.2406.04655}.
#' @references Vehtari A, Gelman A, Gabry J (2017). "Practical Bayesian Model
#' Evaluation Using Leave-One-out Cross-Validation and WAIC."
#' *Statistics and Computing*, **27**(5), 1413-1432. ISSN 0960-3174.
#' \doi{10.1007/s11222-016-9696-4}.
#' @examples
#' data("sim_stvcPoisson")
#' dat <- sim_stvcPoisson[1:100, ]
#'
#' # Fit a spatial-temporal varying coefficient Poisson GLM
#' mod1 <- stvcGLMexact(y ~ x1 + (x1), data = dat, family = "poisson",
#'                      sp_coords = as.matrix(dat[, c("s1", "s2")]),
#'                      time_coords = as.matrix(dat[, "t_coords"]),
#'                      cor.fn = "gneiting-decay",
#'                      process.type = "multivariate",
#'                      sptParams = list(phi_s = 1, phi_t = 1),
#'                      verbose = FALSE, n.samples = 100)
#' @export
stvcGLMexact <- function(formula, data = parent.frame(), family,
                         sp_coords, time_coords, cor.fn,
                         process.type, sptParams,
                         priors, boundary = 0.5, n.samples,
                         loopd = FALSE, loopd.method = "exact", CV.K = 10,
                         loopd.nMC = 500, verbose = TRUE, ...){

  ##### check for unused args #####
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if (!i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

  ##### family #####
  if(missing(family)){
    stop("Family not specified")
  }else{
    if(!is.character(family)){
      stop("Family must be a character string. Choose from c('poisson',
           'binary', 'binomial').")
    }
    family <- tolower(family)
    if(!family %in% c('poisson', 'binary', 'binomial')){
      stop("Invalid family. Choose from c('poisson', 'binary', 'binomial').")
    }
  }

  ##### formula #####
  if(missing(formula)){
    stop("Formula must be specified")
  }
  if(inherits(formula, "formula")){
    holder <- parseFormula2(formula, data)
    if(family == "binomial"){
      if(!dim(holder[[1L]])[2] == 2){
        stop("Response must be of the form cbind(y, n_trials).")
      }
      y <- as.numeric(holder[[1L]][, 1])
      n.binom <- as.numeric(holder[[1L]][, 2])
    }else{
      y <-  as.numeric(holder[[1L]])
    }
    X <- as.matrix(holder[[2L]])
    X.names <- holder[[3L]]
    X_tilde <- as.matrix(holder[[5L]])
    X_tilde.names <- holder[[6L]]
    if(is.null(X_tilde)){
      stop("Formula does not indicate varying coefficients terms")
    }
  }else{
    stop("Formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(X_tilde)

  if(family == "binary"){n.binom <- rep(1.0, n)}
  if(family == "poisson"){n.binom <- rep(0.0, n)}

  if(family == "poisson"){
    if(any(y < 0)){
      stop("family = 'poisson' but y contains negative values")
    }
    if(any(!floor(y) == y)){
      warning("family = 'poisson' but response contains non-integer values")
    }
  }else if(family == "binary"){
    if(any(y != 0 & y != 1)){
      stop("family = 'binary', response can only be either 0 or 1")
    }
  }else if(family == "binomial"){
    if(any(!floor(n.binom) == n.binom)){
      warning("family = 'binomial' but n_trials contains non-integer values")
    }
    if(any(y > n.binom)){
      stop("Number of successes exceeds number of trials in the data.")
    }
  }

  ## storage mode
  storage.mode(y) <- "double"
  storage.mode(n.binom) <- "double"
  storage.mode(X) <- "double"
  storage.mode(X_tilde) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(r) <- "integer"

  ##### coords #####
  if(!is.matrix(sp_coords)){
    stop("sp_coords must n-by-2 matrix of xy-coordinate locations")
  }
  if(ncol(sp_coords) != 2 || nrow(sp_coords) != n){
    stop("either the sp_coords have more than two columns or, number of rows is
         different than data used in the model formula")
  }

  if(!is.matrix(time_coords)){
    stop("time_coords must n-by-1 matrix of temporal coordinates")
  }
  if(ncol(time_coords) != 1 || nrow(time_coords) != n){
    stop("either the time_coords have more than one column or, number of rows is
         different than data used in the model formula")
  }

  storage.mode(sp_coords) <- "double"
  storage.mode(time_coords) <- "double"

  ##### correlation function #####
  if(missing(cor.fn)){
    stop("cor.fn must be specified")
  }
  if(!cor.fn %in% c("gneiting-decay")){
    stop("cor.fn = '", cor.fn, "' is not a valid option; choose from
         c('gneiting-decay')")
  }

  ##### priors #####
  nu.beta <- 0
  nu.z <- 0
  sigmaSq.xi <- 0
  IW.scale <- 0
  missing.flag <- 0

  if(missing(priors)){
    V.beta <- diag(rep(100.0, p))
    nu.beta <- 2.1
    nu.z <- 2.1
    sigmaSq.xi <- 0.1
    IW.scale <- diag(rep(1.0, r))
  }else{
    names(priors) <- tolower(names(priors))
    if(!'v.beta' %in% names(priors)){
      missing.flag <- missing.flag + 1
      V.beta <- diag(rep(100.0, p))
    }else{
      V.beta <- priors[["v.beta"]]
      if(!is.numeric(V.beta) || length(V.beta) != p^2){
        stop(paste("priors[['V.beta']] must be a ", p, "x", p,
                   " covariance matrix.", sep = ""))
      }
    }
    if(!'nu.beta' %in% names(priors)){
      missing.flag <- missing.flag + 1
      nu.beta <- 2.1
    }else{
      nu.beta <- priors[['nu.beta']]
      if(!is.numeric(nu.beta) || length(nu.beta) != 1){
        stop("priors[['nu.beta']] must be a single numeric value.")
      }
      if(nu.beta < 2.1){
        message("Supplied nu.beta is less than 2.1. Setting it to defaults.")
        nu.beta <- 2.1
      }
    }
    if(!'nu.z' %in% names(priors)){
      missing.flag <- missing.flag + 1
      nu.z <- 2.1
    }else{
      nu.z <- priors[['nu.z']]
      if(!is.numeric(nu.z) || length(nu.z) != 1){
        stop("priors[['nu.z']] must be a single numeric value.")
      }
      if(process.type == 'multivariate'){
        if(nu.z < 1){
          message("Supplied nu.z is less than 0.1. Setting it to defaults.")
          nu.z <- 1
        }
      }else{
        if(nu.z < 2.1){
          message("Supplied nu.z is less than 2.1. Setting it to defaults.")
          nu.z <- 2.1
        }
      }
    }
    if(!'sigmasq.xi' %in% names(priors)){
      missing.flag <- missing.flag + 1
      sigmaSq.xi <- 0.1
    }else{
      sigmaSq.xi <- priors[['sigmasq.xi']]
      if(!is.numeric(sigmaSq.xi) || length(sigmaSq.xi) != 1){
        stop("priors[['sigmasq.xi']] must be a single numeric value.")
      }
      if(sigmaSq.xi < 0){
        stop("priors[['sigmaSq.xi']] must be positive real number.")
      }
    }
    if(process.type == 'multivariate'){
      if(!'iw.scale' %in% names(priors)){
        missing.flag <- missing.flag + 1
        IW.scale <- diag(rep(1.0, r))
      }else{
        IW.scale <- priors[['iw.scale']]
        if(!is.numeric(IW.scale) || length(IW.scale) != r^2){
          stop(paste("priors[['IW.scale']] must be a ", r, "x", r,
                     " covariance matrix.", sep = ""))
        }
      }
    }
    if(missing.flag > 0){
      message("Some priors were not supplied. Using defaults.")
    }
  }

  ## storage mode
  storage.mode(nu.beta) <- "double"
  storage.mode(nu.z) <- "double"
  storage.mode(sigmaSq.xi) <- "double"

  ##### spatial process parameters: gneiting-decay #####

  if(missing(sptParams)){
    stop("sptParams (spatial-temporal process parameters) must be supplied.")
  }
  names(sptParams) <- tolower(names(sptParams))

  if(cor.fn == "gneiting-decay"){
    if(process.type %in% c("independent.shared", "multivariate")){
        phi_s <- 0
        phi_t <- 0
        if(!"phi_s" %in% names(sptParams)){
            stop("phi_s (spatial decay parameter) must be supplied.")
        }
        phi_s <- sptParams[["phi_s"]]
        if(!is.numeric(phi_s) || length(phi_s) != 1){
            stop(paste("phi_s must be a scalar when process.type = ",
            process.type))
        }
        if(phi_s <= 0){
            stop("phi_s (spatial decay) must be a positive real.")
        }
        if(!"phi_t" %in% names(sptParams)){
            stop("phi_t (temporal decay parameter) must be supplied.")
        }
        phi_t <- sptParams[["phi_t"]]
        if(!is.numeric(phi_t) || length(phi_t) != 1){
            stop(paste("phi_t must be a scalar when sptShared =", process.type))
        }
        if(phi_t <= 0){
            stop("phi_t (temporal decay) must be a positive real.")
        }
    }else{
        phi_s <- 0
        phi_t <- 0
        if(!"phi_s" %in% names(sptParams)){
            stop("phi_s (spatial decay parameters) must be supplied.")
        }
        phi_s <- sptParams[["phi_s"]]
        if(!is.numeric(phi_s) || length(phi_s) != r){
            stop("When process.type = 'independent', phi_s must be of length ",
            r)
        }
        if(any(phi_s <= 0)){
            stop("phi_s (spatial decay) must be positive reals.")
        }
        if(!"phi_t" %in% names(sptParams)){
            stop("phi_t (temporal decay parameter) must be supplied.")
        }
        phi_t <- sptParams[["phi_t"]]
        if(!is.numeric(phi_t) || length(phi_t) != r){
            stop("When process.type = 'independent', phi_t must be of length ",
            r)
        }
        if(any(phi_t <= 0)){
            stop("phi_t (temporal decay) must be positive reals.")
        }
    }
  }

  ##### Boundary adjustment parameter #####
  epsilon <- 0

  if(!is.numeric(boundary) || length(boundary) != 1){
    stop("boundary must be a numeric scalar between 0 and 1.")
  }else{
    epsilon <- boundary
    if(epsilon <= 0 || epsilon >= 1){
      message("boundary must be in the interval (0, 1). Using default of 0.5.")
      epsilon <- 0.5
    }
  }
  if(family == "binary"){
    if(epsilon < 0.4){
      message("family = binomial'; boundary < 0.4. Setting boundary = 0.4.")
      epsilon <- 0.4
    }
  }

  ## storage.mode
  storage.mode(epsilon) <- "double"

  ##### Leave-one-out setup #####

  if(loopd){
    if(missing(loopd.method)){
      stop("loopd.method must be specified")
    }else{
      loopd.method <- tolower(loopd.method)
    }
    if(!loopd.method %in% c("exact", "cv")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose from c('exact', 'CV').")
    }
    if(loopd.method == "exact"){
      CV.K <- as.integer(0)
    }
    if(loopd.method == "cv"){
      if(n < 100){
        message("Sample size too low for CV. Finding exact LOO-PD.")
        loopd.method <- "exact"
        CV.K <- as.integer(0)
      }else{
        if(CV.K < 10){
          message("CV.K must be at least 10. Setting it to 10.")
          CV.K <- 10
        }else if(CV.K > 20){
          message("CV.K must be at most 20. Setting it to 20.")
          CV.K <- 20
        }
        if(floor(CV.K) != CV.K){
          message("CV.K must be integer. Setting it to nearest integer.")
        }
      }
    }
    if(loopd.nMC < 500){
      message("Number of Monte Carlo samples too low. Using defaults.")
      loopd.nMC = 500
    }
  }else{
    loopd.method <- "none"
  }

  storage.mode(CV.K) <- "integer"
  storage.mode(loopd.nMC) <- "integer"

  ##### sampling setup #####

  if(missing(n.samples)){
    stop("n.samples must be specified.")
  }

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ##### main function call #####
  ptm <- proc.time()

  if(loopd){
    samps <- .Call("stvcGLMexactLOO", y, X, X_tilde, n, p, r, family, n.binom,
                   sp_coords, time_coords, cor.fn,
                   V.beta, nu.beta, nu.z, sigmaSq.xi, IW.scale,
                   process.type, phi_s, phi_t, epsilon,
                   n.samples, loopd, loopd.method, CV.K, loopd.nMC, verbose)
  }else{
    samps <- .Call("stvcGLMexact", y, X, X_tilde, n, p, r, family, n.binom,
                   sp_coords, time_coords, cor.fn,
                   V.beta, nu.beta, nu.z, sigmaSq.xi, IW.scale,
                   process.type, phi_s, phi_t, epsilon,
                   n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- y
  out$X <- X
  out$X.names <- X.names
  out$X.stvc.names <- X_tilde.names
  out$family <- family
  out$sp_coords <- sp_coords
  out$time_coords <- time_coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$spt.params <- sptParams
  out$n.samples <- n.samples
  out$priors <- list(mu.beta = rep(0, p), V.beta = V.beta, nu.beta = nu.beta,
                     nu.z = nu.z, sigmaSq.xi = sigmaSq.xi, IW.scale = IW.scale,
                     boundary = epsilon)
  out$samples <- samps[c("beta", "z", "xi")]
  if(loopd){
    out$loopd <- samps[["loopd"]]
  }
  if(cor.fn == 'gneiting-decay'){
    out$model.params <- list(phi_s = phi_s, phi_t = phi_t)
  }
  out$run.time <- run.time

  class(out) <- "stvcGLMexact"

  return(out)


}