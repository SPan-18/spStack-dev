#' Univariate spatial linear mixed model
#'
#' @param formula a symbolic description of the regression model to be fit. See example below.
#' @param data an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spLM_stack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used to model the spatial dependence structure among the observations. Supported covariance model key words are: \code{"exponential"} and \code{"matern"}. See below for details.
#' @param priors a list with each tag corresponding to a parameter name and containing prior details.
#' @param spParams fixed value of spatial process parameters.
#' @param noise_sp_ratio noise-to-spatial variance ratio
#' @param n.samples number of posterior samples to be generated
#' @param verbose logical.
#' @param ... currently no additional argument
#' @export
spLMexact <- function(formula, data = parent.frame(), coords,
                      cor.fn, priors, spParams, noise_sp_ratio,
                      n.samples, verbose = TRUE, ...){

  ####################################################
  ## check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'", i, "' is not an argument")
  }

  ####################################################
  ## formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified!")}

  if(inherits(formula, "formula")){
    holder <- parseFormula(formula, data)
    y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    X.names <- holder[[3]]
  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)

  ## storage mode
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ####################
  ## coords
  ####################
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or,
    number of rows is different than data used in the model formula")
  }

  coords.D <- 0
  coords.D <- iDist(coords)

  ####################################################
  ## correlation function
  ####################################################
  if(missing(cor.fn)){stop("error: cor.fn must be specified")}
  if(!cor.fn %in% c("exponential", "matern")){stop("cor.fn = '", cor.fn, "' is not a valid option; choose from c('exponential', 'matern').")}

  ####################################################
  ## priors
  ####################################################
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.IG <- 0

  if(missing(priors)){

    warning("prior list not supplied, using defaults.")

  }else{

    names(priors) <- tolower(names(priors))

    ## Setup prior for beta
    if("beta.norm" %in% names(priors)){
      beta.Norm <- priors[["beta.norm"]]
      if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
      if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ", p, ".", sep=""))}
      if(length(beta.Norm[[2]]) != p^2 ){stop(paste("error: beta.Norm[[2]] must be a ", p,"x", p," correlation matrix.", sep=""))}
      beta.prior <- "normal"
    }

    ## Setup prior for sigma.sq
    if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
    sigma.sq.IG <- priors[["sigma.sq.ig"]]

    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
    if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}

  }

  ## storage mode
  storage.mode(sigma.sq.IG) <- "double"

  ####################################################
  ## spatial process parameters
  ####################################################
  phi <- 0
  nu <- 0

  if(missing(spParams)){stop("spParams (spatial process parameters) must be supplied.")}

  names(spParams) <- tolower(names(spParams))

  if(!"phi" %in% names(spParams)){stop("phi must be supplied.")}
  phi <- spParams[["phi"]]

  if(!is.numeric(phi) || length(phi) != 1){stop("phi must be a numeric scalar.")}
  if(phi <= 0){stop("phi (decay parameter) must be a positive real number.")}

  if(cor.fn == "matern"){

    if(!"nu" %in% names(spParams)){stop("nu (smoothness parameter) must be supplied.")}
    nu <- spParams[["nu"]]

    if(!is.numeric(nu) || length(nu) != 1){stop("nu must be a numeric scalar.")}
    if(nu <= 0){stop("nu (smoothness parameter) must be a positive real number.")}

  }

  ## storage mode
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"

  ####################################################
  ## noise-to-spatial variance ratio
  ####################################################
  deltasq <- 0

  if(missing(noise_sp_ratio)){
    warning("noise_sp_ratio not supplied. Using noise_sp_ratio = 1.")
    deltasq = 1
  }else{
    deltasq <- noise_sp_ratio
    if(!is.numeric(deltasq) || length(deltasq) != 1){stop("noise_sp_ratio must be a numeric scalar.")}
    if(deltasq <= 0){stop("noise_sp_ratio must be a positive real number.")}
  }

  ## storage mode
  storage.mode(deltasq) <- "double"

  ####################################################
  ## sampling and setup
  ####################################################

  if(missing(n.samples)){stop("n.samples must be specified.")}

  storage.mode(n.samples) <- "integer"
  storage.mode(verbose) <- "integer"

  ####################################################
  ## main function call
  ####################################################

  out <- .Call("spLMexact", y, X, p, n, coords.D,
               beta.prior, beta.Norm, sigma.sq.IG,
               phi, nu, deltasq, cor.fn, n.samples, verbose)


  return(out)

}
