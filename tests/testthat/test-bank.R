test_that("the bank concentrates on candidates the donors support", {
  set.seed(1)
  LL <- cbind(rnorm(40, 0), rnorm(40, -8), rnorm(40, -8))
  b <- fit_bank(LL, lambda = 0.1)
  expect_s3_class(b, "mbmm_bank")
  expect_true(b$converged)
  expect_equal(sum(b$weights), 1, tolerance = 1e-8)
  expect_gt(b$weights[1], 0.8)
})

test_that("the entropy penalty spreads weight", {
  set.seed(2)
  LL <- cbind(rnorm(30, 0), rnorm(30, -3))
  ent <- function(w) -sum(w * log(pmax(w, 1e-300)))
  expect_gt(ent(fit_bank(LL, lambda = 5)$weights),
            ent(fit_bank(LL, lambda = 0)$weights))
})

test_that("a bank needs at least two candidates", {
  expect_error(fit_bank(matrix(0, 5, 1)), "two candidates")
})

test_that("leave-one-out never scores a participant against its own bank", {
  set.seed(3)
  LL <- matrix(rnorm(12 * 4), 12, 4, dimnames = list(paste0("s", 1:12), NULL))
  loo <- bank_loo(LL, lambda = 1)
  expect_equal(dim(loo$weights), c(12L, 4L))
  expect_length(loo$logLik, 12L)
  full <- fit_bank(LL, lambda = 1)$weights
  expect_false(isTRUE(all.equal(loo$weights[1, ], full)))
})

test_that("data-dependent screening is repeated inside each leave-one-out fold", {
  LL <- rbind(
    held_out = c(100, 0, 0, 0),
    d1 = c(0, 10, 0, 0),
    d2 = c(0, 9, 1, 0),
    d3 = c(0, 8, 0, 1))
  loo <- bank_loo(LL, lambda = 1, nu = 1)

  # Candidate 1 is selected only by the held-out participant. It must not enter
  # that participant's donor-screened support.
  expect_equal(unname(loo$weights["held_out", 1]), 0)
  expect_equal(unname(loo$support_size["held_out"]), 1L)
  expect_equal(sum(loo$weights["held_out", ]), 1, tolerance = 1e-8)
})
