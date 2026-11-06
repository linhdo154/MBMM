test_that("a detected measurement scores highest at the trajectory", {
  ll <- obs_loglik(y = rep(1000, 3), mu = c(500, 1000, 2000), sigma = 0.3)
  expect_equal(which.max(ll), 2L)
})

test_that("a censored measurement scores higher the lower the trajectory", {
  ll <- obs_loglik(y = NA, mu = c(10, 50, 200), sigma = 0.5,
                   lod = 100, censored = TRUE)
  expect_true(all(diff(ll) < 0))
})

test_that("censoring needs a limit of detection", {
  expect_error(obs_loglik(NA, 10, 0.5, censored = TRUE), "lod")
})

test_that("class_loglik integrates over the bank", {
  LL <- matrix(c(log(0.2), log(0.8), log(0.5), log(0.5)), 2, 2, byrow = TRUE)
  w  <- c(0.25, 0.75)
  expect_equal(unname(class_loglik(LL, w)),
               log(c(0.25 * 0.2 + 0.75 * 0.8, 0.25 * 0.5 + 0.75 * 0.5)))
})

test_that("arguments recycle to a common length", {
  y <- c(100, 200, 400)
  expect_length(obs_loglik(y, mu = 150, sigma = 0.3), 3L)
  expect_equal(obs_loglik(y, mu = 150, sigma = 0.3),
               obs_loglik(y, mu = rep(150, 3), sigma = rep(0.3, 3)))
  expect_false(anyNA(obs_loglik(y, mu = 150, sigma = 0.3)))
  mixed <- obs_loglik(y, mu = 150, sigma = 0.3, lod = 120,
                      censored = c(FALSE, TRUE, FALSE))
  expect_length(mixed, 3L)
  expect_false(anyNA(mixed))
})
