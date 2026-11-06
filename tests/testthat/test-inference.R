test_that("the likelihood ratio is non-negative with the right df", {
  set.seed(6)
  n <- 100; X <- cbind(1, rnorm(n), rnorm(n))
  l1 <- rnorm(n, 0.5); l0 <- rnorm(n)
  tt <- gate_test(l1, l0, X)
  expect_s3_class(tt, "mbmm_test")
  expect_gte(tt$statistic, -1e-6)
  expect_equal(tt$df, 2L)
  expect_true(tt$p_chisq >= 0 && tt$p_chisq <= 1)
})

test_that("the permutation p-value is valid and bounded below by 1/(B+1)", {
  set.seed(7)
  n <- 60; x <- rnorm(n); X <- cbind(1, x)
  l1 <- rnorm(n, 1.5 * x); l0 <- rnorm(n)
  tt <- gate_test(l1, l0, X, B = 99, seed = 1)
  expect_length(tt$null, 99L)
  expect_true(tt$p_perm >= 1 / 100 && tt$p_perm <= 1)
})

test_that("T_PO is reported without a chi-squared reference", {
  set.seed(8)
  n <- 50; X <- cbind(1, rnorm(n))
  l1 <- rnorm(n, 1); l0 <- rnorm(n)
  tt <- gate_test(l1, l0, X, firth = TRUE, B = 49, seed = 2)
  expect_true(tt$firth)
  expect_true(is.na(tt$p_chisq))
  expect_false(is.na(tt$p_perm))
})

test_that("a null covariate is not rejected", {
  set.seed(9)
  n <- 80; X <- cbind(1, rnorm(n))
  l1 <- rnorm(n); l0 <- rnorm(n)
  expect_gt(gate_test(l1, l0, X, B = 199, seed = 3)$p_perm, 0.05)
})
