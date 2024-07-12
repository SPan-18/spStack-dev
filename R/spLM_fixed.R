#' Univariate spatial linear mixed model
#'
#' @param formula a symbolic description of the regression model to be fit. See example below.
#' @param data an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spLM_stack} is called.
#' @param coords an \eqn{n \times 2}{n x 2} matrix of the observation coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
#' @param cor.fn a quoted keyword that specifies the correlation function used to model the spatial dependence structure among the observations. Supported covariance model key words are: \code{"exponential"} and \code{"matern"}. See below for details.
#' @param priors a list with each tag corresponding to a parameter name and containing prior details.
#' @param spParams fixed value of spatial process parameters.
#' @param n.samples number of posterior samples to be generated
#' @param verbose logical.
#' @param ... currently no additional arguments.
#' @export
spLM_fixed <- function(formula, data = parent.frame(), coords, cor.fn,
                       priors, spParams, n.samples, verbose = TRUE, ...){

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

  ####################################################
  ## correlation function
  ####################################################
  if(missing(cor.fn)){stop("error: cor.fn must be specified")}
  if(!cor.fn %in% c("exponential", "matern")){stop("error: specified cor.fn '", cor.fn, "' is not a valid option;\n
        choose from c('exponential', 'matern').")}

  ####################################################
  ## priors
  ####################################################
  beta.norm <- 0
  beta.prior <- "flat"
  sigma.sq.IG <- 0
  tau.sq.IG <- 0
  nu.Unif <- 0
  phi.Unif <- 0
  nugget <- FALSE







  ####################################################
  ## spatial process parameters
  ####################################################
  if(missing(spParams)){stop("error: spatial process parameters must be specified")}

  return(list(y, X, X.names))

}
