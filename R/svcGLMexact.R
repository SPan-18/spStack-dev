#' Bayesian spatially varying coefficients generalized linear model
#'
#' @description Fits a Bayesian generalized linear model with spatially varying
#' coefficients under fixed values of the spatial process parameters and some
#' auxiliary model parameters. The output contains posterior samples of the
#' fixed effects, the spatially varying coefficients and, if required, finds
#' leave-one-out predictive densities.
#' @details With this function, we fit a Bayesian hierarchical spatially
#' varying coefficients generalized linear model by sampling exactly from the
#' joint posterior distribution utilizing the generalized conjugate
#' multivariate distribution theory (Bradley and Clinch 2024). Suppose
#' \eqn{\chi = (s_1, \ldots, s_n)} denotes the \eqn{n} spatial locations in
#' \eqn{\mathcal{S} \subset \mathbb{R}^2} at which the response \eqn{y} is
#' observed, with \eqn{y(s)} endowed with a probability law from the natural
#' exponential family
#' \deqn{
#' y(s) \sim \mathrm{EF}(x(s)^\top \beta + \tilde{x}(s)^\top z(s); b(s), \psi)
#' }
#' for some positive parameter \eqn{b(s) > 0} and unit log partition function
#' \eqn{\psi}. Here, \eqn{\tilde{x}(s)} denotes the \eqn{r} covariates with
#' spatially varying coefficients. The supported families (`'poisson'`,
#' `'binomial'` and `'binary'`) and the hierarchical model are those of
#' [stvcGLMexact()], with
#' \deqn{
#' z_j \mid \sigma^2_{z_j} \sim N(0, \sigma^2_{z_j} R(\chi; \phi_j, \nu_j)),
#' \quad \sigma^2_{z_j} \sim \mathrm{IG}(\nu_z/2, \nu_z/2), \quad j = 1,
#' \ldots, r,
#' }
#' where \eqn{R(\chi; \phi, \nu)} is the exponential or Matern spatial
#' correlation matrix with decay \eqn{\phi} and smoothness \eqn{\nu}. If
#' `process.type = 'independent'`, each process has its own \eqn{(\phi_j,
#' \nu_j)}; if `process.type = 'independent.shared'`, they share one. If
#' `process.type = 'multivariate'`, \eqn{z = (z_1, \ldots, z_r)} is a
#' multivariate Gaussian process with correlation matrix \eqn{R(\chi; \phi,
#' \nu)} and an \eqn{r \times r}{r x r} cross-covariance matrix with an
#' inverse-Wishart prior. We fix the process parameters and the
#' hyperparameters \eqn{V_\beta}, \eqn{\nu_\beta}, \eqn{\nu_z} and
#' \eqn{\sigma^2_\xi}; their defaults are as in [stvcGLMexact()].
#' @param formula a symbolic description of the regression model to be fit.
#' Variables in parenthesis are assigned spatially varying coefficients. See
#' examples.
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{svcGLMexact} is called.
#' @param family Specifies the distribution of the response as a member of the
#' exponential family. Supported options are `'poisson'`, `'binomial'` and
#' `'binary'`.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation
#' coordinates in \eqn{\mathbb{R}^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used
#' to model the spatial dependence structure among the observations. Supported
#' covariance model key words are: \code{'exponential'} and \code{'matern'}.
#' @param process.type a quoted keyword specifying the model for the spatial
#' processes: `'independent'`, `'independent.shared'` or `'multivariate'`. See
#' details.
#' @param spParams fixed values of the spatial process parameters, a list with
#' tags `phi` (spatial decay) and, if `cor.fn = 'matern'`, `nu` (smoothness).
#' If `process.type = 'independent'`, each is a vector of length \eqn{r},
#' otherwise a scalar.
#' @param priors (optional) a list with each tag corresponding to a
#' hyperparameter name and containing hyperprior details. Valid tags include
#' `V.beta`, `nu.beta`, `nu.z`, `sigmaSq.xi` and `IW.scale`. Values of `nu.beta`
#' and `nu.z` must be at least 2.1 (`nu.z`, the inverse-Wishart degrees of
#' freedom if `process.type = 'multivariate'`, at least 1). If not supplied,
#' uses defaults.
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
#' @return An object of class \code{svcGLMexact}, which is a list with the
#'  following tags -
#' \describe{
#' \item{priors}{details of the priors used, containing the values of the
#' boundary adjustment parameter (`boundary`), the variance parameter of the
#' fine-scale variation term (`sigmaSq.xi`) and others.}
#' \item{samples}{a list of length 3, containing posterior samples of fixed
#'  effects (\code{beta}), the spatially varying coefficients (\code{z}, an
#'  \eqn{nr \times} \code{n.samples} matrix whose rows \eqn{(j-1)n + 1, \ldots,
#'  jn} correspond to the \eqn{j}-th varying coefficient) and the fine-scale
#'  variation term (\code{xi}).}
#' \item{loopd}{If \code{loopd=TRUE}, contains leave-one-out predictive
#' densities.}
#' \item{model.params}{Values of the fixed parameters: \code{phi} and
#' \code{nu} (\code{NA} for the exponential correlation function).}
#' \item{diagnostics}{a list of diagnostics. Element \code{numerical} is a
#' data frame with one row (one row per process, if \code{process.type =
#' 'independent'}) and columns \code{min.pivot} (the smallest relative
#' Cholesky pivot of the \eqn{n \times n}{n x n} spatial correlation matrix;
#' values below 1e-8 indicate a nearly singular correlation matrix),
#' \code{min.cor} and \code{max.cor} (the correlations of the two
#' farthest-apart and of the two closest locations; values of \code{min.cor}
#' above 0.95 suggest an effective range far exceeding the extent of the data,
#' values of \code{max.cor} below 0.05 nearly uncorrelated locations). If
#' \code{verbose = TRUE}, a "Diagnostics" section is printed when any threshold
#' is crossed.}
#' }
#' The return object might include additional data that can be used for
#' subsequent prediction and/or model fit evaluation.
#' @seealso [svcGLMstack()], [svcLMexact()], [spGLMexact()], [stvcGLMexact()]
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @references Bradley JR, Clinch M (2024). "Generating Independent Replicates
#' Directly from the Posterior Distribution for a Class of Spatial Hierarchical
#' Models." *Journal of Computational and Graphical Statistics*, **0**(0), 1-17.
#' \doi{10.1080/10618600.2024.2365728}.
#' @references Pan S, Zhang L, Bradley JR, Banerjee S (2025). "Bayesian
#' Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking."
#' *Bayesian Analysis*, **In Press**. \doi{10.1214/25-BA1582}.
#' @references Vehtari A, Gelman A, Gabry J (2017). "Practical Bayesian Model
#' Evaluation Using Leave-One-out Cross-Validation and WAIC."
#' *Statistics and Computing*, **27**(5), 1413-1432. ISSN 0960-3174.
#' \doi{10.1007/s11222-016-9696-4}.
#' @examples
#' data(simSpatial)
#' dat <- simSpatial[1:150, ]
#'
#' mod1 <- svcGLMexact(y_svc_pois ~ x1 + x2 + (x1), data = dat,
#'                     family = "poisson",
#'                     coords = as.matrix(dat[, c("s1", "s2")]),
#'                     cor.fn = "matern", process.type = "independent",
#'                     spParams = list(phi = c(8, 4), nu = c(1.5, 1.5)),
#'                     n.samples = 100, verbose = FALSE)
#'
#' # rows 1, ..., n of z hold the varying intercept, rows n+1, ..., 2n the
#' # varying slope of x1
#' n <- nrow(dat)
#' z_hat <- apply(mod1$samples$z, 1, median)
#' cor(z_hat[n + 1:n], dat$z_true)
#' @export
svcGLMexact <- function(formula, data = parent.frame(), family, coords,
                        cor.fn, process.type, spParams, priors,
                        boundary = 0.5, n.samples, loopd = FALSE,
                        loopd.method = "exact", CV.K = 10, loopd.nMC = 500,
                        verbose = TRUE, ...){

  ##### check for unused args #####
  check_dots(...)

  ##### family, formula, process type #####
  dd <- vc_glm_data(formula, data, if(missing(family)) stop("Family not specified") else family)
  family <- dd$family
  n <- dd$n
  p <- dd$p
  r <- dd$r
  process.type <- check_vc_glm_process_type(process.type)
  nR <- vc_nR(process.type, r)

  ##### coords #####
  if(missing(coords)){
    stop("coords must be supplied.")
  }
  if(is.data.frame(coords)){
    coords <- as.matrix(coords)
  }
  if(!is.matrix(coords) || ncol(coords) != 2 || nrow(coords) != n){
    stop("coords must be an n x 2 matrix of spatial coordinates, with n = ", n,
         " the number of observations in the model formula.")
  }
  storage.mode(coords) <- "double"
  check_no_missing(y = dd$y, X = dd$X, X_tilde = dd$X_tilde, n.binom = dd$n.binom, coords = coords)
  check_distinct_coords(coords, what = "spatial locations",
                        hint = "Average the observations at a common location, or model them as spatial-temporal data with stvcGLMexact().")

  ##### correlation function and spatial process parameters #####
  cor.fn <- check_svc_cor_fn(if(missing(cor.fn)) NULL else cor.fn)
  sp <- svc_sp_params(if(missing(spParams)) NULL else spParams, cor.fn, nR, process.type)

  ##### priors #####
  pr <- parse_vc_glm_priors(if(missing(priors)) NULL else priors, p, r, process.type)

  ##### Boundary adjustment parameter #####
  if(!is.numeric(boundary) || length(boundary) != 1){
    stop("boundary must be a numeric scalar between 0 and 1.")
  }
  epsilon <- boundary
  if(epsilon <= 0 || epsilon >= 1){
    message("boundary must be in the interval (0, 1). Using default of 0.5.")
    epsilon <- 0.5
  }
  if(epsilon < 0.1){
    message("boundary < 0.1: the latent pseudo-data of observations at the edge of the support (y = 0, or y = trials) are very diffuse (standard deviation about 1/boundary on the linear predictor scale). A larger boundary is recommended.")
  }
  storage.mode(epsilon) <- "double"

  ##### Leave-one-out setup #####
  if(loopd){
    loopd.method <- tolower(loopd.method)
    if(!loopd.method %in% c("exact", "cv")){
      stop("loopd.method = '", loopd.method, "' is not a valid option; choose from c('exact', 'CV').")
    }
    if(loopd.method == "exact"){
      CV.K <- 0L
    }else{
      if(n < 100){
        message("Sample size too low for CV. Finding exact LOO-PD.")
        loopd.method <- "exact"
        CV.K <- 0L
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
          CV.K <- round(CV.K)
        }
      }
    }
    if(loopd.nMC < 500){
      message("Number of Monte Carlo samples too low. Using defaults.")
      loopd.nMC <- 500
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

  # the routines of stvcGLMexact() with a spatial correlation function: no
  # temporal coordinates, and nu in place of phi_t
  phi <- sp$phi
  nu <- svc_c_nu(sp$nu)
  if(loopd){
    CV.update <- as.integer(resolve_CV_update("auto") == "update")
    samps <- .Call(C_stvcGLMexactLOO, dd$y, dd$X, dd$X_tilde, n, p, r, family, dd$n.binom,
                   coords, NULL, cor.fn,
                   pr$V.beta, pr$nu.beta, pr$nu.z, pr$sigmaSq.xi, pr$IW.scale,
                   process.type, phi, nu, epsilon,
                   n.samples, as.integer(loopd), loopd.method, CV.K, loopd.nMC, CV.update,
                   verbose)
  }else{
    samps <- .Call(C_stvcGLMexact, dd$y, dd$X, dd$X_tilde, n, p, r, family, dd$n.binom,
                   coords, NULL, cor.fn,
                   pr$V.beta, pr$nu.beta, pr$nu.z, pr$sigmaSq.xi, pr$IW.scale,
                   process.type, phi, nu, epsilon,
                   n.samples, verbose)
  }

  run.time <- proc.time() - ptm

  out <- list()
  out$y <- dd$y
  out$X <- dd$X
  out$X.names <- dd$X.names
  out$X.svc.names <- dd$X_tilde.names
  out$family <- family
  out$coords <- coords
  out$cor.fn <- cor.fn
  out$process.type <- process.type
  out$n.samples <- n.samples
  out$priors <- list(mu.beta = rep(0, p), V.beta = pr$V.beta, nu.beta = pr$nu.beta,
                     nu.z = pr$nu.z, sigmaSq.xi = pr$sigmaSq.xi, IW.scale = pr$IW.scale,
                     boundary = epsilon)
  out$samples <- samps[c("beta", "z", "xi")]
  if(loopd){
    out$loopd <- samps[["loopd"]]
  }
  out$model.params <- list(phi = sp$phi, nu = sp$nu)
  out$diagnostics <- list(numerical = collect_diagnostics(list(samps)))
  out$run.time <- run.time

  if(verbose){
    print_diagnostics(out$diagnostics,
                      pivot.hint = "nearly coincident locations, a very small decay parameter phi, or a large smoothness nu")
  }

  class(out) <- "svcGLMexact"

  return(out)

}
