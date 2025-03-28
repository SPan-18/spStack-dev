#' Prediction of latent process at new spatial or temporal locations
#'
#' @description A function to sample from the posterior predictive distribution
#' of the latent spatial or spatial-temporal process.
#' @param mod_out an object returned by [stvcGLMexact()], [stvcGLMstack()], etc.
#' @param coords_new a list of new spatial or spatial-temporal coordinates at
#' which the latent process is to be predicted.
#' @param covars_new a list of new covariates at the new spatial or
#' spatial-temporal coordinates. Only required if the model is either a
#' spatially or a spatial-temporally varying coefficient model.
#' @return A list with the following components:
#' \item{samples}{a list of length `nsamples` containing the posterior
#' predictive samples of the latent process at the new locations or times.}
#' \item{coords_new}{a list of new spatial or spatial-temporal coordinates at
#' which the latent process is to be predicted.}
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [stvcGLMexact()], [stvcGLMstack()]
#' @export
posteriorPredict <- function(mod_out, coords_new, covars_new){

    if(missing(mod_out)){
        stop("Model output is missing.")
    }

    if(inherits(mod_out, 'stvcGLMexact')){

        n <- dim(mod_out$X)[1L]                     # number of observations
        p <- length(mod_out$X.names)                # number of fixed effects
        r <- length(mod_out$X.stvc.names)           # number of varying coefficients
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        storage.mode(r) <- "integer"

        family <- mod_out$family

        if(missing(covars_new)){
            stop("Covariates at new locations or times are missing.")
        }
        if(!is.list(covars_new)){
            stop("covars_new must be a list.")
        }
        if(length(covars_new) != 2 && c("fixed", "vc") %in% names(covars_new)){
            stop("covars_new must be a list with tags 'fixed' and 'vc'.")
        }

        X_new <- covars_new[['fixed']]
        X.tilde_new <- covars_new[['vc']]
        if(!is.matrix(X_new) || !is.matrix(X.tilde_new)){
            stop("Both X_new and X.tilde_new must be matrices.")
        }
        if(nrow(X_new) != nrow(X.tilde_new)){
            stop("X_new and X.tilde_new must have the same number of rows.")
        }
        if(ncol(X_new) != p){
            stop("Number of columns in covariates with fixed effects does not
                 match the number of fixed effects in the model.")
        }
        if(ncol(X.tilde_new) != r){
            stop("Number of columns in columns in covariates with varying
                 coefficients does not match the number of varying coefficients
                 in the model.")
        }
        storage.mode(X_new) <- "double"
        storage.mode(X.tilde_new) <- "double"

        n_pred <- nrow(X_new)                       # number of new predictions
        storage.mode(n_pred) <- "integer"

        # Read spatial-temporal coordinates
        sp_coords <- mod_out$sp_coords
        time_coords <- mod_out$time_coords
        storage.mode(sp_coords) <- "double"
        storage.mode(time_coords) <- "double"

        if(missing(coords_new)){
            stop("Spatial-temporal coordinates at new locations or times are
                 missing.")
        }
        if(!is.list(coords_new)){
            stop("coords_new must be a list.")
        }
        if(length(coords_new) != 2 && c("sp", "time") %in% names(coords_new)){
            stop("coords_new must be a list with tags 'sp' and 'time'.")
        }
        sp_coords_new <- coords_new[['sp']]
        time_coords_new <- as.matrix(coords_new[['time']])
        if(nrow(sp_coords_new) != nrow(time_coords_new)){
            stop("Dimension mismatch in number of new spatial and temporal
                 coordinates.")
        }
        if(ncol(sp_coords_new) != 2){
            stop("Spatial coordinates must have two columns.")
        }
        if(ncol(time_coords_new) != 1){
            stop("Time coordinates must be a vector or a matrix with one
                 column.")
        }
        storage.mode(sp_coords_new) <- "double"
        storage.mode(time_coords_new) <- "double"

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"

        cor.fn <- mod_out$cor.fn
        process.type <- mod_out$process.type
        sptParams <- mod_out$spt.params
        if(cor.fn == "gneiting-decay"){
            phi_s <- sptParams[["phi_s"]]
            phi_t <- sptParams[["phi_t"]]
        }

        # Read samples
        beta_samps <- mod_out$samples[['beta']]
        z_samps <- mod_out$samples[['z']]
        if(all(c('sigmasq.beta', 'z.scale') %in% names(mod_out$samples))){
            sigmasq.beta_samps <- mod_out$samples[['sigmasq.beta']]
            z.scale_samps <- mod_out$samples[['z.scale']]
        }else{
            post_scale <- recoverGLMscale(mod_out)
            sigmasq.beta_samps <- post_scale[['sigmasq.beta']]
            z.scale_samps <- post_scale[['z.scale']]
        }

        # storage mode
        storage.mode(beta_samps) <- "double"
        storage.mode(z_samps) <- "double"
        storage.mode(sigmasq.beta_samps) <- "double"
        storage.mode(z.scale_samps) <- "double"

        # Call C++ function
        samples <- .Call("predict_stvcGLMexact", n, n_pred, p, r, family,
                         X_new, X.tilde_new,
                         sp_coords, time_coords, sp_coords_new, time_coords_new,
                         process.type, cor.fn, phi_s, phi_t, n.samples,
                         beta_samps, z_samps, sigmasq.beta_samps, z.scale_samps)

  }
}