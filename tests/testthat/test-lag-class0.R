test_that("lag integration marginalises and beats a wrong fixed clock", {
  decay <- function(t, th) th$peak * exp(-th$rate * pmax(t, 0))
  cand  <- data.frame(peak = 1000, rate = 0.3)
  shift <- 3
  tt <- c(0, 2, 5, 9)
  obs <- data.frame(id = 1, time = tt,
                    value = decay(tt + shift, cand))
  fixed  <- trajectory_loglik(obs, cand, decay, sigma = 0.2)
  lagged <- trajectory_loglik(obs, cand, decay, sigma = 0.2, lags = 0:5)
  expect_gt(lagged[1, 1], fixed[1, 1])
})

test_that("a single zero lag reproduces no lag at all", {
  decay <- function(t, th) th$peak * exp(-th$rate * t)
  cand <- data.frame(peak = c(100, 400), rate = c(.2, .4))
  obs  <- data.frame(id = rep(1:2, each = 3), time = rep(c(1, 3, 6), 2),
                     value = c(90, 50, 20, 380, 200, 60))
  expect_equal(trajectory_loglik(obs, cand, decay, 0.3),
               trajectory_loglik(obs, cand, decay, 0.3, lags = 0))
})

test_that("the lag prior is honoured", {
  decay <- function(t, th) th$peak * exp(-th$rate * t)
  cand <- data.frame(peak = 1000, rate = 0.3)
  obs  <- data.frame(id = 1, time = c(0, 3), value = decay(c(0, 3), cand))
  flat  <- trajectory_loglik(obs, cand, decay, 0.3, lags = c(0, 5))
  tilt  <- trajectory_loglik(obs, cand, decay, 0.3, lags = c(0, 5),
                             lag_prior = c(0.99, 0.01))
  expect_gt(tilt[1, 1], flat[1, 1])
})

test_that("every censored measurement contributes the same background constant", {
  obs <- data.frame(id = c(1, 2), time = c(1, 1), value = NA_real_,
                    lod = c(100, 10000), censored = TRUE)
  b <- background_loglik(obs, sigma = 1.19)
  expect_equal(unname(b[1]), unname(b[2]))
})

test_that("background_loglik needs a limit of detection", {
  obs <- data.frame(id = 1, time = 1, value = 10, censored = FALSE)
  expect_error(background_loglik(obs, sigma = 1), "lod")
})
