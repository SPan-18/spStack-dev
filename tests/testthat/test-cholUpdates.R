test_that("cholUpdateRankOne works", {
  n <- 10
  A <- matrix(rnorm(n^2), n, n)
  A <- crossprod(A)
  cholA <- chol(A)
  v <- 1:n
  APlusvvT <- A + tcrossprod(v)
  cholA1 <- t(chol(APlusvvT))
  cholA2 <- cholUpdateRankOne(cholA, v, lower = F)
  expect_equal(cholA1, cholA2)
})

test_that("cholUpdateDel works", {
  n <- 10
  A <- matrix(rnorm(n^2), n, n)
  A <- crossprod(A)
  cholA <- chol(A)
  ind <- 2
  A1 <- A[-ind, -ind]
  cholA1 <- t(chol(A1))
  cholA2 <- cholUpdateDel(cholA, del.index = ind, lower = F)
  expect_equal(cholA1, cholA2)
})

test_that("cholUpdateDelBlock works", {
  n <- 10
  A <- matrix(rnorm(n^2), n, n)
  A <- crossprod(A)
  cholA <- chol(A)
  start_ind <- 2
  end_ind <- 6
  del_ind <- c(start_ind:end_ind)
  A1 <- A[-del_ind, -del_ind]
  cholA1 <- t(chol(A1))
  cholA2 <- cholUpdateDelBlock(cholA, start_ind, end_ind, lower = F)
  expect_equal(cholA1, cholA2)
})