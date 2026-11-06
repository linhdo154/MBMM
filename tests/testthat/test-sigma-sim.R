test_that("profiling recovers the scale the data were generated at", {
  set.seed(2)
  decay <- function(t, th) th$peak * exp(-th$rate * t)
  tt <- c(1, 3, 6, 10)
  truth <- data.frame(peak = 1000, rate = 0.25)
  obs <- do.call(rbind, lapply(1:25, function(i) {
    mu <- decay(tt, truth)
    data.frame(id = i, time = tt,
               value = 10^(log10(mu) + rnorm(length(tt), 0, 0.30)))
  }))
  cand <- expand.grid(peak = c(600, 1000, 1600), rate = c(.15, .25, .35))
  p <- profile_sigma(obs, cand, decay, sigma_grid = c(0.1, 0.2, 0.3, 0.5, 0.9))
  expect_true(p$sigma %in% c(0.2, 0.3))
  expect_s3_class(p$profile, "data.frame")
  expect_equal(nrow(p$profile), 5L)
})

test_that("the simulation separates the latent state from any detection", {
  skip_if_not_installed("deSolve")
  s <- simulate_mbmm(n = 40, beta_ab = -1, p_abortive = 0.5,
                     n_pool = 10, seed = 3)
  expect_equal(nrow(s$truth), 40L)
  expect_true(all(c("censored", "lod") %in% names(s$obs)))
  detected <- tapply(!s$obs$censored, s$obs$id, any)
  expect_true(any(detected & s$truth$z == 0))   # contamination exists
  expect_true(all(s$truth$abortive[s$truth$z == 1] == FALSE))
})

test_that("the simulated prevalence tracks the target", {
  skip_if_not_installed("deSolve")
  s <- simulate_mbmm(n = 400, beta_ab = -1, prevalence = 0.3,
                     n_pool = 6, seed = 4)
  expect_lt(abs(mean(s$truth$z) - 0.3), 0.08)
})
