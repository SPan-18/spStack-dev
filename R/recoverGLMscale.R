#' Recover posterior samples of scale parameters of spatial/spatial-temporal
#' generalized linear models
#'
#' @description A function to recover posterior samples of scale parameters that
#' were marginalized out during model fit. This is only applicable for spatial
#' or, spatial-temporal generalized linear models. This function applies on
#' outputs of functions that fits a spatial/spatial-temporal generalized linear
#' model, such as [spGLMexact()], [spGLMstack()], [stvcGLMexact()], and
#' [stvcGLMstack()].
#' @param mod_out an object returned by a fitting a spatial or spatial-temporal
#' GLM.
#' @return An object of the same class as input, and updates the list tagged
#' `samples` with the posterior samples of the scale parameters. The new tags
#' are `sigmasq.beta` and `z.scale`.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [spGLMexact()], [spGLMstack()], [stvcGLMexact()], [stvcGLMstack()]
#' @examples
#' set.seed(1234)
#' data("simPoisson")
#' dat <- simPoisson[1:100, ]
#' mod1 <- spGLMstack(y ~ x1, data = dat, family = "poisson",
#'                    coords = as.matrix(dat[, c("s1", "s2")]), cor.fn = "matern",
#'                    params.list = list(phi = c(3, 5, 7), nu = c(0.5, 1.5),
#'                                       boundary = c(0.5)),
#'                    n.samples = 100,
#'                    loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
#'                    verbose = TRUE)
#'
#' # Recover posterior samples of scale parameters
#' mod1.1 <- recoverGLMscale(mod1)
#'
#' # sample from the stacked posterior distribution
#' post_samps <- stackedSampler(mod1.1)
#' @export
recoverGLMscale <- function(mod_out){

  if(inherits(mod_out, 'stvcGLMexact')){

    n <- dim(mod_out$X)[1L]                     # number of observations
    p <- length(mod_out$X.names)                # number of fixed effects
    r <- length(mod_out$X.stvc.names)           # number of varying coefficients
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(r) <- "integer"

    # Read spatial-temporal coordinates
    sp_coords <- mod_out$sp_coords
    time_coords <- mod_out$time_coords
    cor.fn <- mod_out$cor.fn
    storage.mode(sp_coords) <- "double"
    storage.mode(time_coords) <- "double"

    # Initialize hyperparameters
    nu.beta <- 0
    nu.z <- 0
    sigmaSq.xi <- 0
    V.beta <- matrix(0, p, p)
    IW.scale <- matrix(0, r, r)

    # Read hyperparameters
    priors <- mod_out$priors
    mu.beta <- priors[['mu.beta']]
    V.beta <- priors[['V.beta']]
    nu.beta <- priors[['nu.beta']]
    nu.z <- priors[['nu.z']]
    IW.scale <- priors[['IW.scale']]
    ## storage mode
    storage.mode(nu.beta) <- "double"
    storage.mode(nu.z) <- "double"
    storage.mode(V.beta) <- "double"
    storage.mode(IW.scale) <- "double"

    process.type <- mod_out$process.type
    sptParams <- mod_out$spt.params
    if(cor.fn == "gneiting-decay"){
      phi_s <- sptParams[["phi_s"]]
      phi_t <- sptParams[["phi_t"]]
    }

    n.samples <- mod_out$n.samples
    storage.mode(n.samples) <- "integer"

    beta_samps <- mod_out$samples[['beta']]
    z_samps <- mod_out$samples[['z']]

    scale_samps <- .Call("recoverScale_stvcGLM", n, p, r,
                         sp_coords, time_coords, cor.fn,
                         mu.beta, V.beta, nu.beta, nu.z, IW.scale,
                         process.type, phi_s, phi_t, n.samples,
                         beta_samps, z_samps)

    mod_out$samples[['sigmasq.beta']] <- scale_samps[['sigmasq.beta']]
    mod_out$samples[['z.scale']] <- scale_samps[['z.scale']]

    return(mod_out)

  }else if(inherits(mod_out, 'stvcGLMstack')){

    n <- dim(mod_out$X)[1L]                     # number of observations
    p <- length(mod_out$X.names)                # number of fixed effects
    r <- length(mod_out$X.stvc.names)           # number of varying coefficients
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(r) <- "integer"

    nModels <- mod_out$n.models

    # Read spatial-temporal coordinates
    sp_coords <- mod_out$sp_coords
    time_coords <- mod_out$time_coords
    cor.fn <- mod_out$cor.fn
    storage.mode(sp_coords) <- "double"
    storage.mode(time_coords) <- "double"

    # Initialize hyperparameters
    nu.beta <- 0
    nu.z <- 0
    sigmaSq.xi <- 0
    V.beta <- matrix(0, p, p)
    IW.scale <- matrix(0, r, r)

    # Read hyperparameters
    priors <- mod_out$priors
    mu.beta <- priors[['mu.beta']]
    V.beta <- priors[['V.beta']]
    nu.beta <- priors[['nu.beta']]
    nu.z <- priors[['nu.z']]
    IW.scale <- priors[['IW.scale']]
    ## storage mode
    storage.mode(nu.beta) <- "double"
    storage.mode(nu.z) <- "double"
    storage.mode(V.beta) <- "double"
    storage.mode(IW.scale) <- "double"

    n.samples <- mod_out$n.samples
    storage.mode(n.samples) <- "integer"

    process.type <- mod_out$process.type

    for(i in 1:nModels){
        if(cor.fn == "gneiting-decay"){
            phi_s <- mod_out$candidate.models[[i]][["phi_s"]]
            phi_t <- mod_out$candidate.models[[i]][["phi_t"]]
        }
        beta_samps <- mod_out$samples[[i]][['beta']]
        z_samps <- mod_out$samples[[i]][['z']]

        scale_samps <- .Call("recoverScale_stvcGLM", n, p, r,
                             sp_coords, time_coords, cor.fn,
                             mu.beta, V.beta, nu.beta, nu.z, IW.scale,
                             process.type, phi_s, phi_t, n.samples,
                             beta_samps, z_samps)

        mod_out$samples[[i]][['sigmasq.beta']] <- scale_samps[['sigmasq.beta']]
        mod_out$samples[[i]][['z.scale']] <- scale_samps[['z.scale']]
    }

    return(mod_out)

  }else if(inherits(mod_out, 'spGLMexact')){

    n <- dim(mod_out$X)[1L]                     # number of observations
    p <- length(mod_out$X.names)                # number of fixed effects
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"

    # Read spatial coordinates
    sp_coords <- mod_out$coords
    cor.fn <- mod_out$cor.fn
    coords.D <- 0
    coords.D <- iDist(sp_coords)
    storage.mode(coords.D) <- "double"

    # Initialize hyperparameters
    nu.beta <- 0
    nu.z <- 0
    V.beta <- matrix(0, p, p)

    # Read hyperparameters
    priors <- mod_out$priors
    mu.beta <- priors[['mu.beta']]
    nu.beta <- priors[['nu.beta']]
    nu.z <- priors[['nu.z']]
    V.beta <- priors[['V.beta']]

    ## storage mode
    storage.mode(mu.beta) <- "double"
    storage.mode(nu.beta) <- "double"
    storage.mode(nu.z) <- "double"
    storage.mode(V.beta) <- "double"

    spParams <- mod_out$model.params
    phi <- 0.0
    nu <- 0.0
    if(cor.fn == "matern"){
      phi <- spParams[['phi']]
      nu <- spParams[['nu']]
    }else if(cor.fn == "exponential"){
      phi <- spParams[['phi']]
    }
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"

    n.samples <- mod_out$n.samples
    storage.mode(n.samples) <- "integer"

    beta_samps <- mod_out$samples[['beta']]
    z_samps <- mod_out$samples[['z']]

    scale_samps <- .Call("recoverScale_spGLM", n, p,
                         coords.D, cor.fn,
                         mu.beta, V.beta, nu.beta, nu.z,
                         phi, nu, n.samples,
                         beta_samps, z_samps)

    mod_out$samples[['sigmasq.beta']] <- scale_samps[['sigmasq.beta']]
    mod_out$samples[['sigmasq.z']] <- scale_samps[['sigmasq.z']]

    return(mod_out)

  }else if(inherits(mod_out, 'spGLMstack')){

    n <- dim(mod_out$X)[1L]                     # number of observations
    p <- length(mod_out$X.names)                # number of fixed effects
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"

    # Read spatial coordinates
    sp_coords <- mod_out$coords
    cor.fn <- mod_out$cor.fn
    coords.D <- 0
    coords.D <- iDist(sp_coords)
    storage.mode(coords.D) <- "double"

    # Initialize hyperparameters
    nu.beta <- 0
    nu.z <- 0
    V.beta <- matrix(0, p, p)

    # Read hyperparameters
    priors <- mod_out$priors
    mu.beta <- priors[['mu.beta']]
    nu.beta <- priors[['nu.beta']]
    nu.z <- priors[['nu.z']]
    V.beta <- priors[['V.beta']]

    ## storage mode
    storage.mode(mu.beta) <- "double"
    storage.mode(nu.beta) <- "double"
    storage.mode(nu.z) <- "double"
    storage.mode(V.beta) <- "double"

    n.samples <- mod_out$n.samples
    storage.mode(n.samples) <- "integer"

    nModels <- mod_out$n.models

    for(i in 1:nModels){

      phi <- 0.0
      nu <- 0.0
      if(cor.fn == "matern"){
        phi <- as.numeric(mod_out$candidate.models[i, 'phi'])
        nu <- as.numeric(mod_out$candidate.models[i, 'nu'])
      }else if(cor.fn == "exponential"){
        phi <- as.numeric(mod_out$candidate.models[i, 'phi'])
      }
      storage.mode(phi) <- "double"
      storage.mode(nu) <- "double"

      beta_samps <- mod_out$samples[[i]][['beta']]
      z_samps <- mod_out$samples[[i]][['z']]

      scale_samps <- .Call("recoverScale_spGLM", n, p,
                           coords.D, cor.fn,
                           mu.beta, V.beta, nu.beta, nu.z,
                           phi, nu, n.samples,
                           beta_samps, z_samps)

      mod_out$samples[[i]][['sigmasq.beta']] <- scale_samps[['sigmasq.beta']]
      mod_out$samples[[i]][['sigmasq.z']] <- scale_samps[['sigmasq.z']]

    }

    return(mod_out)

  }

}