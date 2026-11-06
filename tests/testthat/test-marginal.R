test_that("pool_marginal is a log mean-exp within each shared point", {
  ll <- matrix(c(-1, -2, -3, -4, -5, -6), nrow = 2)   # 2 x 3
  out <- pool_marginal(ll, c("a", "a", "b"))
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(out[1, "a"], log(mean(exp(c(-1, -3)))), ignore_attr = TRUE)
  expect_equal(out[2, "a"], log(mean(exp(c(-2, -4)))), ignore_attr = TRUE)
  expect_equal(out[, "b"], c(-5, -6), ignore_attr = TRUE)
})

test_that("a singleton pool leaves the value alone", {
  ll <- matrix(c(-1.5, -2.5), nrow = 1)
  expect_equal(as.vector(pool_marginal(ll, c("x", "y"))), c(-1.5, -2.5))
})

test_that("pool_marginal averages, so it cannot exceed the maximum", {
  set.seed(1)
  ll <- matrix(rnorm(40), nrow = 4)
  sh <- rep(1:2, each = 5)
  out <- pool_marginal(ll, sh)
  mx  <- cbind(apply(ll[, 1:5], 1, max), apply(ll[, 6:10], 1, max))
  expect_true(all(out <= mx + 1e-12))
})

test_that("explicit equal weights reproduce the default", {
  set.seed(2)
  ll <- matrix(rnorm(24), nrow = 3)
  sh <- rep(1:4, each = 2)
  expect_equal(pool_marginal(ll, sh),
               pool_marginal(ll, sh, weights = rep(1, 8)))
})

test_that("pool weights are renormalised within a shared point", {
  ll <- matrix(c(-1, -3), nrow = 1)
  out <- pool_marginal(ll, c("a", "a"), weights = c(3, 1))
  expect_equal(as.vector(out), log(0.75 * exp(-1) + 0.25 * exp(-3)))
})

test_that("column order follows the factor levels", {
  ll <- matrix(c(-1, -2, -3), nrow = 1)
  f  <- factor(c("b", "a", "b"), levels = c("b", "a"))
  expect_equal(colnames(pool_marginal(ll, f)), c("b", "a"))
})

test_that("pool_marginal rejects a mismatched shared vector", {
  expect_error(pool_marginal(matrix(0, 2, 3), c("a", "b")), "ncol")
})

test_that("a zero offset returns the kept columns unchanged", {
  E <- matrix(log(1:8), nrow = 1)
  expect_equal(as.vector(jitter_marginal(E, n_grid = 4, keep = 2:3, offsets = 0L)),
               log(c(2, 3, 6, 7)))
})

test_that("jitter_marginal averages on the natural scale", {
  E <- matrix(log(c(1, 2, 3, 4, 5)), nrow = 1)
  out <- jitter_marginal(E, n_grid = 5, keep = 2:4, offsets = c(-1, 0, 1))
  expect_equal(as.vector(out),
               log(c(mean(1:3), mean(2:4), mean(3:5))))
})

test_that("jitter_marginal keeps the grid coordinate running fastest", {
  E <- matrix(log(1:8), nrow = 1)                      # n_grid 4, K 2
  out <- jitter_marginal(E, n_grid = 4, keep = 2:3, offsets = c(-1, 1))
  expect_equal(as.vector(out),
               log(c(mean(c(1, 3)), mean(c(2, 4)),
                     mean(c(5, 7)), mean(c(6, 8)))))
})

test_that("jitter weights are honoured", {
  E <- matrix(log(c(1, 2, 3)), nrow = 1)
  out <- jitter_marginal(E, n_grid = 3, keep = 2L, offsets = c(-1, 1),
                         weights = c(3, 1))
  expect_equal(as.vector(out), log(0.75 * 1 + 0.25 * 3))
})

test_that("leaving the grid is an error, not a silent wrap", {
  E <- matrix(0, 1, 5)
  expect_error(jitter_marginal(E, 5, keep = 1:5, offsets = c(-1, 1)),
               "leaves the grid")
})

test_that("jitter_marginal needs n_grid to divide the columns", {
  expect_error(jitter_marginal(matrix(0, 1, 7), 4, 1:2, 0L))
})

test_that("the two-level build is identifiable where selection is not", {
  # A flat participant.  Class 1's pool mostly rises; class 0's is mostly flat.
  # Both pools contain an exactly flat member, so SELECTION ties -- each class
  # has a candidate that fits perfectly.  INTEGRATION does not tie, because it
  # asks how a typical member of the pool does, and class 1's typical member
  # predicts a rise the participant did not have.
  tt  <- c(0, 1, 2, 3)
  obs <- data.frame(id = "p", time = tt, value = rep(100, 4))
  rise <- 10^seq(0, 1.5, length.out = 20)      # 1 .. 31.6, mostly rising
  flat <- 10^seq(0, 0.04, length.out = 20)     # 1 .. 1.1, essentially flat
  cand <- rbind(data.frame(gain = rise, g = 1),
                data.frame(gain = flat, g = 0))
  traj <- function(t, th) 100 * th$gain^(t / max(tt))
  LL  <- trajectory_loglik(obs, cand, traj, sigma = 0.1)

  sel <- tapply(LL[1, ], cand$g, max)
  int <- as.vector(pool_marginal(LL, cand$g))
  expect_lt(abs(diff(sel)), 1e-8)              # selection cannot tell them apart
  expect_gt(int[1] - int[2], 1)                # integration prefers class 0
})
