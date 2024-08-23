#' Sample from the *stacked posterior* distribution
#'
#' @description A helper function to sample from the stacked posterior
#' distribution to obtain final posterior samples for analysis. This function
#' applies on outputs of functions like [spLMstack()].
#' @param mod_out output object from [spLMstack()].
#' @param n.samples number of posterior samples to draw from the stacked
#' posterior. If it exceeds the number of posterior draws used in the original
#' function, then it throws a warning and resamples within the available
#' samples. It is recommended, to run the original function with the number of
#' samples that is desired at the final step. If missing, inherits the number
#' of posterior samples from the orinial output.
#' @return An object of class \code{stacked_posterior}, which is a list with the
#' following tags -
#' \item{beta}{samples of the fixed effect from the stacked joint posterior.}
#' \item{sigmaSq}{samples of the variance parameter from the stacked joint
#' posterior.}
#' \item{z}{samples of the spatial random effects from the stacked joint
#' posterior.}
#' @details After obtaining the optimal stacking weights
#' \eqn{\hat{w}_1, \ldots, \hat{w}_G}, posterior inference of quantities of
#' interest subsequently proceed from the *stacked posterior*,
#' \deqn{
#' \tilde{p}(\cdot \mid y) = \sum_{g = 1}^G \hat{w}_g p(\cdot \mid y, M_g),
#' }
#' where \eqn{\mathcal{M} = \{M_1, \ldots, M_g\}} is the collection of candidate
#' models.
#' @seealso [spLMstack()]
#' @export
stackedSampler <- function(mod_out, n.samples){

    if(inherits(mod_out, 'spLMstack')){

      n_obs <- dim(mod_out$X)[1L]
      p_obs <- dim(mod_out$X)[2L]
      n_post <- dim(mod_out$samples[[1L]][['beta']])[2]

      if(missing(n.samples)){
        n.samples <- n_post
        ids <- sample(1:n_post, size = n.samples, replace = FALSE)
      }else{
        if(n.samples > n_post){
        warning("Number of samples required exceeds number of posterior samples.
                To prevent resampling, run spLMstack() with higher n.samples.")
        ids <- sample(1:n_post, size = n.samples, replace = TRUE)
        }else{
        ids <- sample(1:n_post, size = n.samples, replace = FALSE)
        }
      }

      post_samples <- sapply(1:n.samples, function(x){
        model_id <- sample(1:mod_out$n.models, 1,
                           prob = mod_out$stacking.weights)
        return(c(mod_out$samples[[model_id]]$beta[, ids[x]],
                 mod_out$samples[[model_id]]$sigmaSq[ids[x]],
                 mod_out$samples[[model_id]]$z[, ids[x]]))
            })
      stacked_samps <- list(beta = post_samples[1:p_obs, ],
                            sigmaSq = post_samples[(p_obs + 1), ],
                            z = post_samples[(p_obs + 1) + 1:n_obs, ])
      rownames(stacked_samps[['beta']]) = mod_out$X.names

    }else{
      stop("Invalid model output class. Input must be an output from
           spLMstack().")
      # Append to this list as new functions are added
    }

    class(stacked_samps) <- "stacked_posterior"
    return(stacked_samps)
}