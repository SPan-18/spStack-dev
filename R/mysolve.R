#' Solve a positive-definite linear system
#'
#' Solves the linear system \eqn{Ax = b}{Ax = b} where \eqn{A}{A} is a positive-definite matrix and \eqn{b}{b} is a vector.
#'
#' @useDynLib spStack
#' @param A an \eqn{n\times n}{nxn} positive-definite matrix
#' @param b an \eqn{n\times 1}{nx1} matrix/vector
#' @returns The \eqn{n\times 1}{nx1} solution to the system
#' @export
#' @examples
#' A <- rbind(c(1, 2, 3),
#'            c(2, 5, 7),
#'            c(3, 7, 14))
#' b <- c(1, 0, 1)
#' (x <- mysolve(A, b)) # 5 -2 0
#'
#' @keywords utilities
mysolve <- function(A, b) {

    n <- nrow(A)
    res <- matrix(0, n)

    storage.mode(A) <- "double"
    storage.mode(b) <- "double"
    storage.mode(res) <- "double"
    storage.mode(n) <- "integer"

    .Call("mysolveC", A, b, n)

}
