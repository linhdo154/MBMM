test_that("mbmm_class_loglik favors the bank whose curves match the data", {
  obs <- data.frame(time = c(0, 1, 2), y = c(1, 2, 3))

  close_bank <- data.frame(id = "donorA", bundle_id = 1, time = c(0, 1, 2),
                            pred = c(1.05, 1.95, 3.1))
  far_bank <- data.frame(id = "donorB", bundle_id = 1, time = c(0, 1, 2),
                          pred = c(5, 5, 5))

  ll_close <- mbmm_class_loglik(obs, close_bank, sigma = 0.2)
  ll_far <- mbmm_class_loglik(obs, far_bank, sigma = 0.2)
  expect_gt(ll_close, ll_far)
})

test_that("mbmm_class_loglik excludes the donor's own particles (leave-one-out)", {
  obs <- data.frame(time = c(0, 1), y = c(1, 1))
  bank <- data.frame(id = c("S1", "S2"), bundle_id = c(1, 2),
                      time = c(0, 0), pred = c(1, 1))
  # only S2's row overlaps t=1 too -> add it
  bank <- rbind(bank, data.frame(id = c("S1", "S2"), bundle_id = c(1, 2),
                                  time = c(1, 1), pred = c(1, 1)))

  ll_excl <- mbmm_class_loglik(obs, bank, sigma = 0.5, exclude_donor = "S1")
  ll_incl <- mbmm_class_loglik(obs, bank, sigma = 0.5, exclude_donor = NULL)
  expect_true(is.finite(ll_excl))
  expect_true(is.finite(ll_incl))
})

test_that("mbmm_class_loglik returns -Inf when there is no time overlap", {
  obs <- data.frame(time = c(100, 200), y = c(1, 1))
  bank <- data.frame(id = "A", bundle_id = 1, time = c(0, 1), pred = c(1, 1))
  expect_equal(mbmm_class_loglik(obs, bank, sigma = 0.3), -Inf)
})

test_that("mbmm_class_loglik supports a target column and per-target sigma", {
  obs <- data.frame(time = c(0, 0), target = c("IgG", "IgM"), y = c(2, 4))
  bank <- data.frame(id = "A", bundle_id = 1, time = c(0, 0),
                      target = c("IgG", "IgM"), pred = c(2, 4))
  sigma <- c(IgG = 0.1, IgM = 0.5)
  ll <- mbmm_class_loglik(obs, bank, sigma = sigma)
  expect_true(is.finite(ll))
})

test_that("mbmm_score_subjects produces one row/column per subject/class", {
  obs_list <- list(
    S1 = data.frame(time = 0:1, y = c(1, 1)),
    S2 = data.frame(time = 0:1, y = c(5, 5))
  )
  bankA <- data.frame(id = "donor", bundle_id = 1, time = 0:1, pred = c(1, 1))
  bankB <- data.frame(id = "donor", bundle_id = 1, time = 0:1, pred = c(5, 5))

  logL <- mbmm_score_subjects(obs_list, list(A = bankA, B = bankB), sigma = 0.3,
                               verbose = FALSE)
  expect_equal(nrow(logL), 2)
  expect_setequal(names(logL), c("id", "logL_A", "logL_B"))
  # S1 should be much better explained by bank A, S2 by bank B
  expect_gt(logL$logL_A[logL$id == "S1"], logL$logL_B[logL$id == "S1"])
  expect_gt(logL$logL_B[logL$id == "S2"], logL$logL_A[logL$id == "S2"])
})
