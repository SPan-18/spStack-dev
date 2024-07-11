test_that("mysolve works", {
  A <- matrix(rnorm(25), ncol = 5)
  A <- crossprod(A)
  b <- 1:5
  expect_equal(mysolve(A, b), as.numeric(solve(A) %*% b))
})
