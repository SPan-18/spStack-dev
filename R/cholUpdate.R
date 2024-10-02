#' Different Cholesky factor updates
#'
#'
#' @description Provides functions that implements different types of updates of
#' a Cholesky factor that includes rank-one update, single row/column deletion
#' update and a block deletion update.
#' @name cholUpdate
#' @param A an \eqn{n\times n} triangular matrix
#' @param v an \eqn{n\times 1} matrix/vector
#' @param alpha scalar; if not supplied, default is 1
#' @param beta scalar; if not supplied, default is 1
#' @param del.index an integer from 1 to \eqn{n} indicating the row/column to be
#' deleted
#' @param del.start an integer from 1 to \eqn{n} indicating the first row/column
#' of a block to be deleted, must be at least 1 less than `del.end`
#' @param del.end an integer from 1 to \eqn{n} indicating the last row/column
#' of a block to be deleted, must be at least 1 more than `del.start`
#' @param lower logical; if `A` is lower-triangular or not
#' @returns An \eqn{m \times m} lower-triangular `matrix`` with \eqn{m = n} in
#' case of `cholUpdateRankOne()`, \eqn{m = n - 1} in case of `cholUpdateDel()`,
#' and, \eqn{m = n - n_k} in case of `cholUpdateDelBlock()` where \eqn{n_k} is
#' the size of the block removed.
#' @details Suppose \eqn{B = AA^\top} is a \eqn{n \times n} matrix with \eqn{A}
#' being its lower-triangular Cholesky factor. Then rank-one update corresponds
#' to finding the Cholesky factor of the matrix
#' \eqn{C = \alpha B + \beta vv^\top} for some \eqn{\alpha,\beta\in\mathbb{R}}
#' given \eqn{A} (see, Krause and Igel 2015). Similarly, single row/column
#' deletion update corresponds to finding the Cholesky factor of the
#' \eqn{(n-1)\times(n-1)} matrix \eqn{B_i} which is obtained by removing the
#' \eqn{i}-th row and column of \eqn{B}, given \eqn{A} for some
#' \eqn{i - 1, \ldots, n}. Lastly, block deletion corresponds to finding the
#' Cholesky factor of the \eqn{(n-n_k)\times(n-n_k)} matrix \eqn{B_{I}} for a
#' subset \eqn{I} of \eqn{\{1, \ldots, n\}} containing \eqn{n_k} consecutive
#' indices, given the factor \eqn{A}.
#' @references Oswin Krause and Christian Igel. 2015. A More Efficient Rank-one
#' Covariance Matrix Update for Evolution Strategies. In *Proceedings of the
#' 2015 ACM Conference on Foundations of Genetic Algorithms XIII* (FOGA '15).
#' Association for Computing Machinery, New York, NY, USA, 129-136.
#' \doi{10.1145/2725494.2725496}.
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @examples
#' n <- 10
#' A <- matrix(rnorm(n^2), n, n)
#' A <- crossprod(A)
#' cholA <- chol(A)
#'
#' ## Rank-1 update
#' v <- 1:n
#' APlusvvT <- A + tcrossprod(v)
#' cholA1 <- t(chol(APlusvvT))
#' cholA2 <- cholUpdateRankOne(cholA, v, lower = FALSE)
#' print(all(abs(cholA1 - cholA2) < 1E-9))
#'
#' ## Single Row-deletion update
#' ind <- 2
#' A1 <- A[-ind, -ind]
#' cholA1 <- t(chol(A1))
#' cholA2 <- cholUpdateDel(cholA, del.index = ind, lower = FALSE)
#' print(all(abs(cholA1 - cholA2) < 1E-9))
#'
#' ## Block-deletion update
#' start_ind <- 2
#' end_ind <- 6
#' del_ind <- c(start_ind:end_ind)
#' A1 <- A[-del_ind, -del_ind]
#' cholA1 <- t(chol(A1))
#' cholA2 <- cholUpdateDelBlock(cholA, start_ind, end_ind, lower = FALSE)
#' print(all(abs(cholA1 - cholA2) < 1E-9))
NULL

#' @rdname cholUpdate
#' @export
cholUpdateRankOne <- function(A, v, alpha, beta, lower = TRUE){

  n <- nrow(A)
  if(length(v) != n){ stop("Dimension mismatch.") }
  if(missing(alpha)){ alpha <- 1 }
  if(missing(beta)){ beta <- 1 }

  storage.mode(A) <- "double"
  storage.mode(v) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(alpha) <- "double"
  storage.mode(beta) <- "double"
  storage.mode(lower) <- "integer"

  .Call("R_cholRankOneUpdate", A, n, v, alpha, beta, lower)

}

#' @rdname cholUpdate
#' @export
cholUpdateDel <- function(A, del.index, lower = TRUE){

  n <- nrow(A)
  if(missing(del.index)){
    stop("del.index missing.")
  }
  if(n < 3){ stop("Must be a matrix of dimension at least 3.")}
  if(del.index < 0 || del.index > n){
    stop("Index to delete out of bounds.")
  }

  storage.mode(A) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(del.index) <- "integer"
  storage.mode(lower) <- "integer"

  .Call("R_cholRowDelUpdate", A, n, del.index, lower)

}

#' @rdname cholUpdate
#' @export
cholUpdateDelBlock <- function(A, del.start, del.end, lower = TRUE){

  n <- nrow(A);
  if(missing(del.start)){ stop("Start index missing.") }
  if(missing(del.end)){ stop("End index missing.") }
  if(del.start > del.end || del.start == del.end){ stop("Start index must be at least 1 less than end index.") }
  if(del.start < 0 || del.end > n){ stop("Indices to delete out of bounds.") }
  n_del <- del.end - del.start + 1
  if(n_del > n - 2){ stop("Maximum block size to delete is ", n - 2, ".") }

  storage.mode(A) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(del.start) <- "integer"
  storage.mode(del.end) <- "integer"
  storage.mode(lower) <- "integer"

  .Call("R_cholRowBlockDelUpdate", A, n, del.start, del.end, lower)

}
