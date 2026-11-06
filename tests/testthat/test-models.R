test_that("amplitudes and fluxes round-trip", {
  f <- amp_to_flux(400, 1200, 50, 3000, 2500, ms = 0.3)
  a <- flux_to_amp(f$A0, f$Ws, f$Wl, f$dWs, f$dWl, f$ms)
  expect_equal(unlist(a[c("A0", "C1", "C2", "D1", "D2")]),
               c(A0 = 400, C1 = 1200, C2 = 50, D1 = 3000, D2 = 2500),
               tolerance = 1e-8)
})

test_that("the antibody trajectory starts at baseline and is continuous at the boost", {
  tb <- 84 / 30
  a <- function(t) andraud_trajectory(t, 400, 1200, 0, 3000, 2500, ms = 0.3)
  expect_equal(a(0), 400, tolerance = 1e-6)
  expect_equal(a(tb - 1e-7), a(tb + 1e-7), tolerance = 1e-3)
})

test_that("the trajectory approaches the post-boost plateau", {
  expect_equal(andraud_trajectory(400, 400, 1200, 0, 3000, 2500, ms = 0.3),
               3000, tolerance = 1)
})

test_that("the class coordinate changes sign with the reference level", {
  expect_gt(andraud_levels(400, 1200, 3000, "C1"), 0)
  expect_gt(andraud_levels(400, 1200, 3000, "A0"), 0)
  expect_lt(andraud_levels(400, 4000, 3000, "C1"), 0)
})

test_that("every eclipse candidate establishes by construction", {
  d <- eclipse_design(200, seed = 1)
  R0 <- d$beta * d$T0 * d$p / (d$delta * d$c_clear)
  expect_true(all(R0 > 1))
  expect_true(all(d$Td > 0))
})
