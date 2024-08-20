#' @importFrom stats is.empty.model model.matrix model.response terms
parseFormula <- function(formula, data, intercept = TRUE, justX = FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept) {
        attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    X <- as.matrix(X)  # X matrix
    xvars <- dimnames(X)[[2L]]  # X variable names
    xobs <- dimnames(X)[[1L]]  # X observation names
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

# internal function: input a list of candidate values of all model parameters
# expands list of vectors, called inside spLMstack and similar functions
candidate_models <- function(params_list){

    models <- expand.grid(params_list)
    models_list <- apply(models, 1, function(x) as.vector(x, mode = "list"))

    return(models_list)

}