#' Sample from the stacked posterior distribution
#'
#' @description A helper function to sample from the stacked posterior
#' distribution to obtain final posterior samples that can be used for
#' subsequent analysis. This function applies on outputs of functions
#' [spLMstack()] and [spGLMstack()].
#' @param mod_out an object that is an output of a model fit or a prediction
#' task, i.e., the class should be either `spLMstack`, 'pp.spLMstack',
#' `spGLMstack`, `pp.spGLMstack`, `stvcGLMexact`, or `pp.stvcGLMexact`.
#' @param n.samples (optional) If missing, inherits the number
#' of posterior samples from the original output. Otherwise, it specifies
#' number of posterior samples to draw from the stacked posterior. If it exceeds
#' the number of posterior draws used in the original function, then a message
#' is thrown and the samples are obtained by resampling. We recommended running
#' the original model fit/prediction with enough samples.
#' @return An object of class \code{stacked_posterior}, which is a list that
#' includes the following tags -
#' \describe{
#' \item{beta}{samples of the fixed effect from the stacked joint posterior.}
#' \item{z}{samples of the spatial random effects from the stacked joint
#' posterior.}
#' }
#' The list may also include other scale parameters corresponding to the model.
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
#' set.seed(1234)
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
              To prevent resampling, run *stack() with higher n.samples.")
      ids <- sample(1:n_post, size = n.samples, replace = TRUE)
    }else{
    ids <- sample(1:n_post, size = n.samples, replace = FALSE)
    }
  }

  if(inherits(mod_out, c('spLMstack', 'pp.spLMstack',
                         'spGLMstack', 'pp.spGLMstack',
                         'stvcGLMstack', 'pp.stvcGLMstack'))){

    nModels <- mod_out$n.models
    model_id <- sample(seq_len(nModels), n.samples, replace = TRUE, prob = mod_out$stacking.weights)

    param_names <- names(mod_out$samples[[1L]])
    result <- vector("list", length(param_names))
    names(result) <- param_names

    for(param in param_names){

      # Determine shape of the parameter
      first_param <- mod_out$samples[[1]][[param]]
      is_matrix <- is.matrix(first_param)
      d <- if (is_matrix) nrow(first_param) else 1

      # Preallocate
      if(is_matrix){
        result[[param]] <- matrix(NA, nrow = d, ncol = n.samples)
      }else{
        result[[param]] <- numeric(n.samples)
      }

      for(i in seq_len(n.samples)){
        m <- model_id[i]
        s <- ids[i]
        val <- mod_out$samples[[m]][[param]]
        if(is.matrix(val)){
          result[[param]][, i] <- val[, s]
        }else{
          result[[param]][i] <- val[s]
        }
      }

    }

    if("beta" %in% names(result)){
      rownames(result[['beta']]) <- mod_out$X.names
    }

  }else{
    stop("Invalid model output class. Input must be an output from either of the
         following functions: spLMstack(), spGLMstack(), stvcGLMstack().")
  }

  class(result) <- "stacked_posterior"
  return(result)

}