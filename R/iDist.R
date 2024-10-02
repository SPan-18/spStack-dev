#' Calculate distance matrix
#' 
#' Computes the inter-site Euclidean distance matrix for one or two sets of
#' points.
#'
#' @usage iDist(coords.1, coords.2, ...)
#' @param coords.1 an \eqn{n\times p} matrix with each row corresponding to
#' a point in \eqn{p}-dimensional space.
#' @param coords.2 an \eqn{m\times p} matrix with each row corresponding to
#' a point in \eqn{p} dimensional space. If this is missing then
#' \code{coords.1} is used.
#' @param ... currently no additional arguments.
#' @returns The \eqn{n\times n} or \eqn{n\times m} inter-site Euclidean
#' distance matrix.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @examples
#' n <- 10
#' p1 <- cbind(runif(n),runif(n))
#' m <- 5
#' p2 <- cbind(runif(m),runif(m))
#' D <- iDist(p1, p2)
#' @keywords utilities
#' @export
iDist <- function(coords.1, coords.2, ...) {

    if (!is.matrix(coords.1))
        coords.1 <- as.matrix(coords.1)

    if (missing(coords.2))
        coords.2 <- coords.1

    if (!is.matrix(coords.2))
        coords.2 <- as.matrix(coords.2)

    if (ncol(coords.1) != ncol(coords.2))
        stop("error: ncol(coords.1) != ncol(coords.2)")

    p <- ncol(coords.1)
    n1 <- nrow(coords.1)
    n2 <- nrow(coords.2)


    D <- matrix(0, n1, n2)

    storage.mode(coords.1) <- "double"
    storage.mode(coords.2) <- "double"
    storage.mode(D) <- "double"
    storage.mode(n1) <- "integer"
    storage.mode(n2) <- "integer"
    storage.mode(p) <- "integer"

    .Call("idist", coords.1, n1, coords.2, n2, p, D)

    return(D)
}
