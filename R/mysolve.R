#' Solves the linear system Ax = b where A is positive definite
#'
#' @useDynLib spStackdev
#' @param A an \eqn{n\times n}{nxn} positive definite matrix
#' @param b an \eqn{n\times 1}{nx1} matrix/vector
#' @returns The \eqn{n\times 1}{nx1} solution to the system
#' @export
#' @keywords utilities
mysolve <- function(A, b){

  n <- nrow(A)
  res <- matrix(0, n)

  storage.mode(A) <- "double"
  storage.mode(b) <- "double"
  storage.mode(res) <- "double"
  storage.mode(n) <- "integer"

  .Call("mysolveC", A, b, n)

}
