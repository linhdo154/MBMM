test_that("restriction renormalises and reports the mass", {
  w <- c(0.1, 0.2, 0.3, 0.4)
  r <- restrict_bank(w, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(r$mass, 0.3)
  expect_equal(sum(r$weights), 1)
  expect_equal(r$weights[1] / r$weights[2], w[1] / w[2])
  expect_equal(r$weights[3:4], c(0, 0))
})

test_that("class prevalence cancels after restriction", {
  w1 <- c(0.4, 0.4, 0.1, 0.1); w2 <- c(0.1, 0.1, 0.4, 0.4)
  reg <- c(TRUE, TRUE, FALSE, FALSE)
  expect_equal(restrict_bank(w1, reg)$weights, restrict_bank(w2, reg)$weights)
})

test_that("an empty region is an error and a thin one warns", {
  expect_error(restrict_bank(c(0.5, 0.5), c(FALSE, FALSE)), "no population mass")
  expect_warning(restrict_bank(c(0.999, 0.001), c(FALSE, TRUE)), "amplifies")
})
