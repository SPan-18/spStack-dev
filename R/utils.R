#' @importFrom stats is.empty.model model.matrix model.response terms
parseFormula <- function(formula, data, intercept = TRUE, justX = FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept) {
        attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    X <- as.matrix(X)  # X matrix
    xvars <- dimnames(X)[[2]]  # X variable names
    xobs <- dimnames(X)[[1]]  # X observation names
    if (justX) {
        Y <- NULL
    } else {
        Y <- as.matrix(model.response(mf, "numeric"))  # Y matrix
    }

    return(list(Y, X, xvars, xobs))

}

# internal function: checks if input is integer
is_integer <- function(x) {
    is.numeric(x) && (floor(x) == x)
}

#' Create candidate models
#'
#' @description `candidate_models` creates an `stack.model.list` object
#' that is passed into the functions `spLMstack` or `spGLMstack`.
#' @details Supply a list with each entry containing candidate values of the
#' parameter of the same name as that entry. See example below.
#' @param params_list a list containing candidate values of parameters
#' required for the model that you want to fit. For names of the parameters,
#' look up the function that implements the model. See example below.
#' @return a `stack.model.list` object
#' @examples
#' mod_list <- list(phi = c(3, 5, 10),
#'                  nu = c(0.5, 1, 1.5))
#' cand_mods <- candidate_models(mod_list)
#' # pass this to spLMstack() in the mod_params argument
#' @export
candidate_models <- function(params_list){

    models <- expand.grid(params_list)
    models_list <- apply(models, 1, function(x) as.vector(x, mode = "list"))

    class(models_list) <- "spStack.model.list"

    return(models_list)
}