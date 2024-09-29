test_that("cholUpdateRankOne", {
  n <- 10
  A <- matrix(rnorm(n^2), n, n)
  A <- crossprod(A)
  cholA <- chol(A)
  v <- 1:n
  APlusvvT <- A + tcrossprod(v)
  cholA1 <- t(chol(APlusvvT))
  cholA2 <- cholUpdateRankOne(cholA, v, lower = F)
  max_diff <- as.numeric(max(abs(cholA1 - cholA2)) < 1E-8)
  expect_equal(max_diff, 1)
})

test_that("cholUpdateDel", {
  n <- 10
  A <- matrix(rnorm(n^2), n, n)
  A <- crossprod(A)
  cholA <- chol(A)
  ind <- 2
  A1 <- A[-ind, -ind]
  cholA1 <- t(chol(A1))
  cholA2 <- cholUpdateDel(cholA, del.index = ind, lower = F)
  max_diff <- as.numeric(max(abs(cholA1 - cholA2)) < 1E-9)
  expect_equal(max_diff, 1)
})

test_that("cholUpdateDelBlock", {
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
  max_diff <- as.numeric(max(abs(cholA1 - cholA2)) < 1E-9)
  expect_equal(max_diff, 1)
})