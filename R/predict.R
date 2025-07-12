#' Prediction of latent process at new spatial or temporal locations
#'
#' @description A function to sample from the posterior predictive distribution
#' of the latent spatial or spatial-temporal process.
#' @param mod_out an object returned by any model fit under fixed
#' hyperparameters or using predictive stacking, i.e., [spLMexact()],
#' [spLMstack()], [spGLMexact()], [spGLMstack()], [stvcGLMexact()], or
#' [stvcGLMstack()].
#' @param coords_new a list of new spatial or spatial-temporal coordinates at
#' which the latent process, the mean, and the response is to be predicted.
#' @param covars_new a list of new covariates at the new spatial or
#' spatial-temporal coordinates. See examples for the structure of this list.
#' @param joint a logical value indicating whether to return the joint posterior
#' predictive samples of the latent process at the new locations or times.
#' Defaults to `FALSE`.
#' @param nBinom_new a vector of the number of trials for each new prediction
#' location or time. Only required if the model family is `"binomial"`. Defaults
#' to a vector of ones, indicating one trial for each new prediction.
#' @return A modified object with the class name preceeded by the identifier
#' `pp` separated by a `.`. For example, if input is of class `spLMstack`, then
#' the output of this prediction function would be `pp.spLMstack`. The entry
#' with the tag `samples` is updated and will include samples from the posterior
#' predictive distribution of the latent process, the mean, and the response at
#' the new locations or times. An entry with the tag `prediction` is added
#' and contains the new coordinates and covariates, and whether the joint
#' posterior predictive samples were requested.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [spLMexact()], [spLMstack()], [spGLMexact()], [spGLMstack()],
#' [stvcGLMexact()], [stvcGLMstack()]
#' @examples
#' set.seed(1234)
#' # training and test data sizes
#' n_train <- 100
#' n_pred <- 10
#'
#' # Example 1: Spatial linear model
#' # load and split data into training and prediction sets
#' data(simGaussian)
#' dat <- simGaussian
#' dat_train <- dat[1:n_train, ]
#' dat_pred <- dat[n_train + 1:n_pred, ]
#'
#' # fit a spatial linear model using predictive stacking
#' mod1 <- spLMstack(y ~ x1, data = dat_train,
#'                   coords = as.matrix(dat_train[, c("s1", "s2")]),
#'                   cor.fn = "matern",
#'                   params.list = list(phi = c(1.5, 3, 5), nu = c(0.75, 1.25),
#'                                      noise_sp_ratio = c(0.5, 1, 2)),
#'                   n.samples = 1000, loopd.method = "psis",
#'                   parallel = FALSE, solver = "ECOS", verbose = TRUE)
#'
#' # prepare new coordinates and covariates for prediction
#' sp_pred <- as.matrix(dat_pred[, c("s1", "s2")])
#' X_new <- as.matrix(cbind(rep(1, n_pred), dat_pred$x1))
#'
#' # carry out posterior prediction
#' mod.pred <- posteriorPredict(mod1, coords_new = sp_pred, covars_new = X_new,
#'                              joint = TRUE)
#'
#' # sample from the stacked posterior and posterior predictive distribution
#' post_samps <- stackedSampler(mod.pred)
#'
#' # analyze posterior samples
#' postpred_z <- post_samps$z.pred
#' post_z_summ <- t(apply(postpred_z, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
#' z_combn <- data.frame(z = dat_pred$z_true, zL = post_z_summ[, 1],
#'                       zM = post_z_summ[, 2], zU = post_z_summ[, 3])
#' library(ggplot2)
#' ggplot(data = z_combn, aes(x = z)) +
#'   geom_errorbar(aes(ymin = zL, ymax = zU), width = 0.05, alpha = 0.15, color = "skyblue") +
#'   geom_point(aes(y = zM), size = 0.25, color = "darkblue", alpha = 0.5) +
#'   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "solid") +
#'   xlab("True z1") + ylab("Posterior of z1") + theme_bw() +
#'   theme(panel.background = element_blank(), aspect.ratio = 1)
#' @export
posteriorPredict <- function(mod_out, coords_new, covars_new, joint = FALSE,
                             nBinom_new){

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
        if(family == "binomial"){
            if(missing(nBinom_new)){
                nBinom_new <- rep(1, n_pred)  # default to 1 trial if not provided
            }
            if(!is.vector(nBinom_new) || length(nBinom_new) != n_pred){
                stop("nBinom_new must be a vector of the same length as the
                     number of new predictions.")
            }
        }else{
            nBinom_new <- rep(0, n_pred)
        }
        storage.mode(nBinom_new) <- "integer"

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

        dup_rows <- find_approx_matches(cbind(sp_coords_new, time_coords_new),
                                        cbind(sp_coords, time_coords))
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatio-temporal coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

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
        if('z.scale' %in% names(mod_out$samples)){
            z.scale_samps <- mod_out$samples[['z.scale']]
        }else{
            mod_out <- recoverGLMscale(mod_out)
            z.scale_samps <- mod_out$samples[['z.scale']]
        }

        # storage mode
        storage.mode(beta_samps) <- "double"
        storage.mode(z_samps) <- "double"
        storage.mode(z.scale_samps) <- "double"
        storage.mode(joint) <- "integer"

        # Call C++ function
        samps <- .Call("predict_stvcGLM", n, n_pred, p, r, family,
                       nBinom_new, X_new, X.tilde_new,
                       sp_coords, time_coords, sp_coords_new, time_coords_new,
                       process.type, cor.fn, phi_s, phi_t, n.samples,
                       beta_samps, z_samps, z.scale_samps,
                       joint)

        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))
        mod_out$samples[['z.pred']] <- samps[['z.pred']]
        mod_out$samples[['mu.pred']] <- samps[['mu.pred']]
        mod_out$samples[['y.pred']] <- samps[['y.pred']]

        class(mod_out) <- "pp.stvcGLMexact"
        return(mod_out)

    }else if(inherits(mod_out, 'stvcGLMstack')){

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
        if(family == "binomial"){
            if(missing(nBinom_new)){
                nBinom_new <- rep(1, n_pred)  # default to 1 trial if not provided
            }
            if(!is.vector(nBinom_new) || length(nBinom_new) != n_pred){
                stop("nBinom_new must be a vector of the same length as the
                     number of new predictions.")
            }
        }else{
            nBinom_new <- rep(0, n_pred)
        }
        storage.mode(nBinom_new) <- "integer"

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

        dup_rows <- find_approx_matches(cbind(sp_coords_new, time_coords_new),
                                        cbind(sp_coords, time_coords))
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatio-temporal coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"
        storage.mode(joint) <- "integer"

        cor.fn <- mod_out$cor.fn
        process.type <- mod_out$process.type

        nModels <- length(mod_out$candidate.models)

        for(i in 1:nModels){

            # read hyperparameters
            if(cor.fn == "gneiting-decay"){
                phi_s <- mod_out$candidate.models[[i]][["phi_s"]]
                phi_t <- mod_out$candidate.models[[i]][["phi_t"]]
            }

            # read samples
            beta_samps <- mod_out$samples[[i]][['beta']]
            z_samps <- mod_out$samples[[i]][['z']]
            if('z.scale' %in% names(mod_out$samples[[i]])){
                z.scale_samps <- mod_out$samples[[i]][['z.scale']]
            }else{
                mod_out <- recoverGLMscale(mod_out)
                z.scale_samps <- mod_out$samples[[i]][['z.scale']]
            }

            # storage mode
            storage.mode(beta_samps) <- "double"
            storage.mode(z_samps) <- "double"
            storage.mode(z.scale_samps) <- "double"

            # Call C++ function
            samps <- .Call("predict_stvcGLM", n, n_pred, p, r, family,
                           nBinom_new, X_new, X.tilde_new,
                           sp_coords, time_coords, sp_coords_new, time_coords_new,
                           process.type, cor.fn, phi_s, phi_t, n.samples,
                           beta_samps, z_samps, z.scale_samps,
                           joint)

            mod_out$samples[[i]][['z.pred']] <- samps[['z.pred']]
            mod_out$samples[[i]][['mu.pred']] <- samps[['mu.pred']]
            mod_out$samples[[i]][['y.pred']] <- samps[['y.pred']]

        }

        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))
        class(mod_out) <- "pp.stvcGLMstack"
        return(mod_out)

    }else if(inherits(mod_out, 'spGLMexact')){

        n <- dim(mod_out$X)[1L]                     # number of observations
        p <- length(mod_out$X.names)                # number of fixed effects
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"

        family <- mod_out$family

        if(missing(covars_new)){
            stop("Covariates at new locations are missing.")
        }
        if(!is.matrix(covars_new)){
            stop("covars_new must be a matrix.")
        }
        if(ncol(covars_new) != p){
            stop("Number of columns in covars_new must match the number of fixed effects in the model.")
        }
        storage.mode(covars_new) <- "double"

        n_pred <- nrow(covars_new)                   # number of new predictions
        storage.mode(n_pred) <- "integer"

        if(family == "binomial"){
            if(missing(nBinom_new)){
                nBinom_new <- rep(1, n_pred)  # default to 1 trial if not provided
            }
            if(!is.vector(nBinom_new) || length(nBinom_new) != n_pred){
                stop("nBinom_new must be a vector of the same length as the
                     number of new predictions.")
            }
        }else{
            nBinom_new <- rep(0, n_pred)
        }
        storage.mode(nBinom_new) <- "integer"

        # Read spatial coordinates
        sp_coords <- mod_out$coords
        storage.mode(sp_coords) <- "double"

        if(missing(coords_new)){
            stop("Spatial coordinates at new locations are missing.")
        }
        if(!is.matrix(coords_new)){
            stop("coords_new must be a matrix.")
        }
        if(nrow(coords_new) != n_pred){
            stop("Number of rows in coords_new must match the number of rows of covars_new.")
        }
        if(ncol(coords_new) != 2){
            stop("Spatial coordinates must have two columns.")
        }
        storage.mode(coords_new) <- "double"

        dup_rows <- find_approx_matches(coords_new, sp_coords)
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatial coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"

        cor.fn <- mod_out$cor.fn

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

        # Read samples
        beta_samps <- mod_out$samples[['beta']]
        z_samps <- mod_out$samples[['z']]
        if('sigmasq.z' %in% names(mod_out$samples)){
            sigmasq.z_samps <- mod_out$samples[['sigmasq.z']]
        }else{
            mod_out <- recoverGLMscale(mod_out)
            sigmasq.z_samps <- mod_out$samples[['sigmasq.z']]
        }

        # storage mode
        storage.mode(beta_samps) <- "double"
        storage.mode(z_samps) <- "double"
        storage.mode(sigmasq.z_samps) <- "double"
        storage.mode(joint) <- "integer"

        # Call C++ function
        samps <- .Call("predict_spGLM", n, n_pred, p, family,
                       nBinom_new, covars_new, sp_coords, coords_new,
                       cor.fn, phi, nu, n.samples,
                       beta_samps, z_samps, sigmasq.z_samps,
                       joint)

        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))
        mod_out$samples[['z.pred']] <- samps[['z.pred']]
        mod_out$samples[['mu.pred']] <- samps[['mu.pred']]
        mod_out$samples[['y.pred']] <- samps[['y.pred']]

        class(mod_out) <- "pp.spGLMexact"
        return(mod_out)

    }else if(inherits(mod_out, 'spGLMstack')){

        n <- dim(mod_out$X)[1L]                     # number of observations
        p <- length(mod_out$X.names)                # number of fixed effects
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"

        family <- mod_out$family

        if(missing(covars_new)){
            stop("Covariates at new locations are missing.")
        }
        if(!is.matrix(covars_new)){
            stop("covars_new must be a matrix.")
        }
        if(ncol(covars_new) != p){
            stop("Number of columns in covars_new must match the number of fixed effects in the model.")
        }
        storage.mode(covars_new) <- "double"

        n_pred <- nrow(covars_new)                   # number of new predictions
        storage.mode(n_pred) <- "integer"

        if(family == "binomial"){
            if(missing(nBinom_new)){
                nBinom_new <- rep(1, n_pred)  # default to 1 trial if not provided
            }
            if(!is.vector(nBinom_new) || length(nBinom_new) != n_pred){
                stop("nBinom_new must be a vector of the same length as the
                     number of new predictions.")
            }
        }else{
            nBinom_new <- rep(0, n_pred)
        }
        storage.mode(nBinom_new) <- "integer"

        # Read spatial coordinates
        sp_coords <- mod_out$coords
        storage.mode(sp_coords) <- "double"

        if(missing(coords_new)){
            stop("Spatial coordinates at new locations are missing.")
        }
        if(!is.matrix(coords_new)){
            stop("coords_new must be a matrix.")
        }
        if(nrow(coords_new) != n_pred){
            stop("Number of rows in coords_new must match the number of rows of covars_new.")
        }
        if(ncol(coords_new) != 2){
            stop("Spatial coordinates must have two columns.")
        }
        storage.mode(coords_new) <- "double"

        dup_rows <- find_approx_matches(coords_new, sp_coords)
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatial coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"
        storage.mode(joint) <- "integer"

        cor.fn <- mod_out$cor.fn

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
            if('sigmasq.z' %in% names(mod_out$samples[[i]])){
                sigmasq.z_samps <- mod_out$samples[[i]][['sigmasq.z']]
            }else{
                mod_out <- recoverGLMscale(mod_out)
                sigmasq.z_samps <- mod_out$samples[[i]][['sigmasq.z']]
            }
            # storage mode
            storage.mode(beta_samps) <- "double"
            storage.mode(z_samps) <- "double"
            storage.mode(sigmasq.z_samps) <- "double"

            # Call C++ function
            samps <- .Call("predict_spGLM", n, n_pred, p,
                           family, nBinom_new, covars_new, sp_coords, coords_new,
                           cor.fn, phi, nu, n.samples,
                           beta_samps, z_samps, sigmasq.z_samps,
                           joint)

            mod_out$samples[[i]][['z.pred']] <- samps[['z.pred']]
            mod_out$samples[[i]][['mu.pred']] <- samps[['mu.pred']]
            mod_out$samples[[i]][['y.pred']] <- samps[['y.pred']]

        }
        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))
        class(mod_out) <- "pp.spGLMstack"
        return(mod_out)

    }else if(inherits(mod_out, 'spLMexact')){

        n <- dim(mod_out$X)[1L]                     # number of observations
        p <- length(mod_out$X.names)                # number of fixed effects
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"

        if(missing(covars_new)){
            stop("Covariates at new locations are missing.")
        }
        if(!is.matrix(covars_new)){
            stop("covars_new must be a matrix.")
        }
        if(ncol(covars_new) != p){
            stop("Number of columns in covars_new must match the number of fixed effects in the model.")
        }
        storage.mode(covars_new) <- "double"

        n_pred <- nrow(covars_new)                   # number of new predictions
        storage.mode(n_pred) <- "integer"

        # Read spatial coordinates
        sp_coords <- mod_out$coords
        storage.mode(sp_coords) <- "double"

        if(missing(coords_new)){
            stop("Spatial coordinates at new locations are missing.")
        }
        if(!is.matrix(coords_new)){
            stop("coords_new must be a matrix.")
        }
        if(nrow(coords_new) != n_pred){
            stop("Number of rows in coords_new must match the number of rows of covars_new.")
        }
        if(ncol(coords_new) != 2){
            stop("Spatial coordinates must have two columns.")
        }
        storage.mode(coords_new) <- "double"

        dup_rows <- find_approx_matches(coords_new, sp_coords)
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatial coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"
        storage.mode(joint) <- "integer"

        cor.fn <- mod_out$cor.fn

        # Read hyperparameters
        phi <- 0.0
        nu <- 0.0
        deltasq <- 0.0
        if(cor.fn == "matern"){
            phi <- mod_out$model.params[['phi']]
            nu <- mod_out$model.params[['nu']]
            deltasq <- mod_out$model.params[['noise_sp_ratio']]
        }else if(cor.fn == "exponential"){
            phi <- mod_out$model.params[['phi']]
            deltasq <- mod_out$model.params[['noise_sp_ratio']]
        }

        # Read samples
        beta_samps <- mod_out$samples[['beta']]
        z_samps <- mod_out$samples[['z']]
        sigmaSq_samps <- mod_out$samples[['sigmaSq']]

        # storage mode
        storage.mode(beta_samps) <- "double"
        storage.mode(z_samps) <- "double"
        storage.mode(sigmaSq_samps) <- "double"
        storage.mode(phi) <- "double"
        storage.mode(nu) <- "double"
        storage.mode(deltasq) <- "double"

        # Call C++ function
        samps <- .Call("predict_spLM", n, n_pred, p,
                       covars_new, sp_coords, coords_new,
                       cor.fn, phi, nu, deltasq,
                       beta_samps, z_samps, sigmaSq_samps,
                       n.samples, joint)

        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))
        mod_out$samples[['z.pred']] <- samps[['z.pred']]
        mod_out$samples[['mu.pred']] <- samps[['mu.pred']]
        mod_out$samples[['y.pred']] <- samps[['y.pred']]

        class(mod_out) <- "pp.spLMexact"
        return(mod_out)

    }else if(inherits(mod_out, 'spLMstack')){

        n <- dim(mod_out$X)[1L]                     # number of observations
        p <- length(mod_out$X.names)                # number of fixed effects
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"

        if(missing(covars_new)){
            stop("Covariates at new locations are missing.")
        }
        if(!is.matrix(covars_new)){
            stop("covars_new must be a matrix.")
        }
        if(ncol(covars_new) != p){
            stop("Number of columns in covars_new must match the number of fixed effects in the model.")
        }
        storage.mode(covars_new) <- "double"

        n_pred <- nrow(covars_new)                   # number of new predictions
        storage.mode(n_pred) <- "integer"

        # Read spatial coordinates
        sp_coords <- mod_out$coords
        storage.mode(sp_coords) <- "double"

        if(missing(coords_new)){
            stop("Spatial coordinates at new locations are missing.")
        }
        if(!is.matrix(coords_new)){
            stop("coords_new must be a matrix.")
        }
        if(nrow(coords_new) != n_pred){
            stop("Number of rows in coords_new must match the number of rows of covars_new.")
        }
        if(ncol(coords_new) != 2){
            stop("Spatial coordinates must have two columns.")
        }
        storage.mode(coords_new) <- "double"

        dup_rows <- find_approx_matches(coords_new, sp_coords)
        if(dup_rows$any_match){
            matched_preview <- dup_rows$matched_rows
            n_show <- min(6, nrow(matched_preview))
            matched_preview <- matched_preview[seq_len(n_show), , drop = FALSE]

            # Convert to character matrix for pretty printing
            matched_str <- apply(matched_preview, 1, function(r) paste(format(r, digits = 6), collapse = ", "))
            matched_str <- paste0("  [", seq_len(n_show), "] ", matched_str, collapse = "\n")

            stop(sprintf(
                "New spatial coordinates match %d existing coordinate(s).\nFirst %d matched row(s):\n%s\nPlease provide new coordinates that do not match existing ones.",
                dup_rows$num_matches, n_show, matched_str))
        }

        n.samples <- mod_out$n.samples
        storage.mode(n.samples) <- "integer"
        storage.mode(joint) <- "integer"

        cor.fn <- mod_out$cor.fn
        nModels <- mod_out$n.models

        for(i in 1:nModels){

            phi <- 0.0
            nu <- 0.0
            deltasq <- 0.0
            if(cor.fn == "matern"){
                phi <- as.numeric(mod_out$candidate.models[i, 'phi'])
                nu <- as.numeric(mod_out$candidate.models[i, 'nu'])
                deltasq <- as.numeric(mod_out$candidate.models[i, 'noise_sp_ratio'])
            }else if(cor.fn == "exponential"){
                phi <- as.numeric(mod_out$candidate.models[i, 'phi'])
                deltasq <- as.numeric(mod_out$candidate.models[i, 'noise_sp_ratio'])
            }
            storage.mode(phi) <- "double"
            storage.mode(nu) <- "double"
            storage.mode(deltasq) <- "double"

            beta_samps <- mod_out$samples[[i]][['beta']]
            z_samps <- mod_out$samples[[i]][['z']]
            sigmaSq_samps <- mod_out$samples[[i]][['sigmaSq']]

            # storage mode
            storage.mode(beta_samps) <- "double"
            storage.mode(z_samps) <- "double"
            storage.mode(sigmaSq_samps) <- "double"

            # Call C++ function
            samps <- .Call("predict_spLM", n, n_pred, p,
                           covars_new, sp_coords, coords_new,
                           cor.fn, phi, nu, deltasq,
                           beta_samps, z_samps, sigmaSq_samps,
                           n.samples, joint)

            mod_out$samples[[i]][['z.pred']] <- samps[['z.pred']]
            mod_out$samples[[i]][['mu.pred']] <- samps[['mu.pred']]
            mod_out$samples[[i]][['y.pred']] <- samps[['y.pred']]

        }

        # Prepare return object
        mod_out$prediction <- list(coords_new = coords_new,
                                   covars_new = covars_new,
                                   joint.pred = as.logical(joint))

        class(mod_out) <- "pp.spLMstack"
        return(mod_out)

    }
    else{
        stop("Input model object must be of class:
        'spLMexact', 'spLMstack', 'spGLMexact', 'spGLMstack', 'stvcGLMexact', 'stvcGLMstack'.")
    }

}