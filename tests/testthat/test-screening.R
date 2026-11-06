test_that("screening keeps each participant's best candidates", {
  set.seed(1)
  LL <- matrix(rnorm(20 * 100), 20, 100)
  s <- screen_candidates(LL, nu = 5)
  expect_lte(s$after, 20 * 5)
  expect_equal(s$before, 100L)
  for (i in 1:20) expect_true(which.max(LL[i, ]) %in% s$keep)
})

test_that("nu beyond the family keeps everything", {
  LL <- matrix(rnorm(10 * 6), 10, 6)
  expect_equal(screen_candidates(LL, nu = 50)$after, 6L)
})

test_that("the positivity screen rejects excursions between observations", {
  cand <- data.frame(a = c(1, 1, 1), b = c(0.1, 2, 0.5))
  f <- function(t, th) th$a - th$b * t
  s <- screen_positive(cand, f, grid = seq(0, 1, by = 0.02), floor = 0.01)
  expect_equal(s$keep, c(TRUE, FALSE, TRUE))
  expect_equal(s$rejected, 1L)
})

test_that("a screen on the observation grid alone misses an excursion", {
  cand <- data.frame(amp = 1, w = 1)
  dip <- function(t, th) ifelse(t > 1.4 & t < 1.6, -1, th$amp)
  expect_true(screen_positive(cand, dip, grid = c(0, 1, 2))$keep)
  expect_false(screen_positive(cand, dip, grid = seq(0, 2, by = 0.05))$keep)
})
