# Different Cholesky factor updates

Provides functions that implements different types of updates of a
Cholesky factor that includes rank-one update, single row/column
deletion update and a block deletion update.

## Usage

``` r
cholUpdateRankOne(A, v, alpha, beta, lower = TRUE)

cholUpdateDel(A, del.index, lower = TRUE)

cholUpdateDelBlock(A, del.start, del.end, lower = TRUE)
```

## Arguments

- A:

  an \\n\times n\\ triangular matrix

- v:

  an \\n\times 1\\ matrix/vector

- alpha:

  scalar; if not supplied, default is 1

- beta:

  scalar; if not supplied, default is 1

- lower:

  logical; if `A` is lower-triangular or not

- del.index:

  an integer from 1 to \\n\\ indicating the row/column to be deleted

- del.start:

  an integer from 1 to \\n\\ indicating the first row/column of a block
  to be deleted, must be at least 1 less than `del.end`

- del.end:

  an integer from 1 to \\n\\ indicating the last row/column of a block
  to be deleted, must be at least 1 more than `del.start`

## Value

An \\m \times m\\ lower-triangular `matrix` with \\m = n\\ in case of
`cholUpdateRankOne()`, \\m = n - 1\\ in case of `cholUpdateDel()`, and,
\\m = n - n_k\\ in case of `cholUpdateDelBlock()` where \\n_k\\ is the
size of the block removed.

## Details

Suppose \\B = AA^\top\\ is a \\n \times n\\ matrix with \\A\\ being its
lower-triangular Cholesky factor. Then rank-one update corresponds to
finding the Cholesky factor of the matrix \\C = \alpha B + \beta
vv^\top\\ for some \\\alpha,\beta\in\mathbb{R}\\ given \\A\\ (see,
Krause and Igel 2015). Similarly, single row/column deletion update
corresponds to finding the Cholesky factor of the \\(n-1)\times(n-1)\\
matrix \\B_i\\ which is obtained by removing the \\i\\-th row and column
of \\B\\, given \\A\\ for some \\i - 1, \ldots, n\\. Lastly, block
deletion corresponds to finding the Cholesky factor of the
\\(n-n_k)\times(n-n_k)\\ matrix \\B\_{I}\\ for a subset \\I\\ of \\\\1,
\ldots, n\\\\ containing \\n_k\\ consecutive indices, given the factor
\\A\\.

## References

Oswin Krause and Christian Igel. 2015. "A More Efficient Rank-one
Covariance Matrix Update for Evolution Strategies". In *Proceedings of
the 2015 ACM Conference on Foundations of Genetic Algorithms XIII* (FOGA
'15). Association for Computing Machinery, New York, NY, USA, 129-136.
[doi:10.1145/2725494.2725496](https://doi.org/10.1145/2725494.2725496) .

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
n <- 10
A <- matrix(rnorm(n^2), n, n)
A <- crossprod(A)
cholA <- chol(A)

## Rank-1 update
v <- 1:n
APlusvvT <- A + tcrossprod(v)
cholA1 <- t(chol(APlusvvT))
cholA2 <- cholUpdateRankOne(cholA, v, lower = FALSE)
print(all(abs(cholA1 - cholA2) < 1E-9))
#> [1] TRUE

## Single Row-deletion update
ind <- 2
A1 <- A[-ind, -ind]
cholA1 <- t(chol(A1))
cholA2 <- cholUpdateDel(cholA, del.index = ind, lower = FALSE)
print(all(abs(cholA1 - cholA2) < 1E-9))
#> [1] TRUE

## Block-deletion update
start_ind <- 2
end_ind <- 6
del_ind <- c(start_ind:end_ind)
A1 <- A[-del_ind, -del_ind]
cholA1 <- t(chol(A1))
cholA2 <- cholUpdateDelBlock(cholA, start_ind, end_ind, lower = FALSE)
print(all(abs(cholA1 - cholA2) < 1E-9))
#> [1] TRUE
```
