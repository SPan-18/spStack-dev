#' Sample from the stacked posterior distribution
#'
#' @description A helper function to sample from the stacked posterior
#' distribution to obtain final posterior samples that can be used for
#' subsequent analysis. This function applies on outputs of functions
#' [spLMstack()] and [spGLMstack()].
#' @param mod_out an object of class `spLMstack` or `spGLMstack`.
#' @param n.samples (optional) If missing, inherits the number
#' of posterior samples from the original output. Otherwise, it specifies
#' number of posterior samples to draw from the stacked posterior. If it exceeds
#' the number of posterior draws used in the original function, then a warning
#' is thrown and the samples are obtained by resampling. It is recommended, to
#' run the original function with enough samples.
#' @return An object of class \code{stacked_posterior}, which is a list that
#' includes the following tags -
#' \describe{
#' \item{beta}{samples of the fixed effect from the stacked joint posterior.}
#' \item{z}{samples of the spatial random effects from the stacked joint
#' posterior.}
#' }
#' In case of model output of class `spLMstack`, the list additionally contains
#' `sigmaSq` which are the samples of the variance parameter from the stacked
#' joint posterior of the spatial linear model. For model output of class
#' `spGLMstack`, the list also contains `xi` which are the samples of the
#' fine-scale variation term from the stacked joint posterior of the spatial
#' generalized linear model.
#' @details After obtaining the optimal stacking weights
#' \eqn{\hat{w}_1, \ldots, \hat{w}_G}, posterior inference of quantities of
#' interest subsequently proceed from the *stacked* posterior,
#' \deqn{
#' \tilde{p}(\cdot \mid y) = \sum_{g = 1}^G \hat{w}_g p(\cdot \mid y, M_g),
#' }
#' where \eqn{\mathcal{M} = \{M_1, \ldots, M_g\}} is the collection of candidate
#' models.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [spLMstack()], [spGLMstack()]
#' @examples
#' data(simGaussian)
#' dat <- simGaussian[1:100, ]
#'
#' mod1 <- spLMstack(y ~ x1, data = dat,
#'                   coords = as.matrix(dat[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   params.list = list(phi = c(1.5, 3),
#'                                      nu = c(0.5, 1),
#'                                      noise_sp_ratio = c(1)),
#'                   n.samples = 1000, loopd.method = "exact",
#'                   parallel = FALSE, solver = "ECOS", verbose = TRUE)
#' print(mod1$solver.status)
#' print(mod1$run.time)
#'
#' post_samps <- stackedSampler(mod1)
#' post_beta <- post_samps$beta
#' print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
#' @export
stackedSampler <- function(mod_out, n.samples){

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

  if(inherits(mod_out, 'spLMstack')){

    post_samples <- sapply(1:n.samples, function(x){
      model_id <- sample(1:mod_out$n.models, 1, prob = mod_out$stacking.weights)
      return(c(mod_out$samples[[model_id]]$beta[, ids[x]],
               mod_out$samples[[model_id]]$sigmaSq[ids[x]],
               mod_out$samples[[model_id]]$z[, ids[x]]))
    })
    stacked_samps <- list(beta = post_samples[1:p_obs, ],
                          sigmaSq = post_samples[(p_obs + 1), ],
                          z = post_samples[(p_obs + 1) + 1:n_obs, ])
    rownames(stacked_samps[['beta']]) = mod_out$X.names

  }else if(inherits(mod_out, 'spGLMstack')){

    post_samples <- sapply(1:n.samples, function(x){
      model_id <- sample(1:mod_out$n.models, 1, prob = mod_out$stacking.weights)
      return(c(mod_out$samples[[model_id]]$beta[, ids[x]],
               mod_out$samples[[model_id]]$z[, ids[x]],
               mod_out$samples[[model_id]]$xi[, ids[x]]))
    })
    stacked_samps <- list(beta = post_samples[1:p_obs, ],
                          z = post_samples[(p_obs + 1:n_obs), ],
                          xi = post_samples[(p_obs + n_obs) + 1:n_obs, ])
    rownames(stacked_samps[['beta']]) = mod_out$X.names

  }else{
    stop("Invalid model output class. Input must be an output from either
         spLMstack() or spGLMstack() functions.")
  }

    class(stacked_samps) <- "stacked_posterior"
    return(stacked_samps)
}