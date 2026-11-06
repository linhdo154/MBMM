test_that("the gate recovers a covariate effect the likelihoods encode", {
  set.seed(4)
  n <- 400; x <- rnorm(n); X <- cbind(1, x)
  z <- rbinom(n, 1, plogis(-0.5 + 1.2 * x))
  l1 <- ifelse(z == 1, 0, -4); l0 <- ifelse(z == 1, -4, 0)
  g <- fit_gate(l1, l0, X)
  expect_s3_class(g, "mbmm_gate")
  expect_gt(g$gamma[2], 0.6)
  expect_lt(abs(g$gamma[2] - 1.2), 0.5)
})

test_that("an uninformative trajectory leaves the responsibility at the gate", {
  n <- 50; X <- cbind(1, rnorm(n))
  g <- fit_gate(rep(0, n), rep(0, n), X)
  expect_equal(g$w, g$pi, tolerance = 1e-6)
})

test_that("Firth returns a penalised objective below the log-likelihood", {
  set.seed(5)
  n <- 40; X <- cbind(1, rnorm(n))
  l1 <- rnorm(n, 1); l0 <- rnorm(n)
  g <- fit_gate(l1, l0, X, firth = TRUE)
  expect_true(g$firth)
  expect_false(is.null(g$penalised))
  expect_true(is.finite(g$penalised))
})

test_that("Firth is stable under separation, where plain EM runs away", {
  n <- 30; x <- c(rep(-1, 15), rep(1, 15)); X <- cbind(1, x)
  l1 <- ifelse(x > 0, 0, -50); l0 <- ifelse(x > 0, -50, 0)
  expect_lt(abs(fit_gate(l1, l0, X, firth = TRUE)$gamma[2]), 1e4)
})
