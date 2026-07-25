# Simple linear trajectory used throughout these tests:
# y(t) = a + b * t
linear_simulate <- function(theta_row, times) {
     as.numeric(theta_row$a) +
          as.numeric(theta_row$b) * times
}

linear_prior <- function(n) {
     data.frame(
          a = stats::runif(n, 0, 5),
          b = stats::runif(n, -2, 2)
     )
}

test_that("mbmm_pool_bank gives every donor equal representation", {
     particles <- data.frame(
          id = rep(c("A", "B", "C"), times = c(10, 5, 8)),
          x = rnorm(23)
     )
     
     # NULL uses the minimum available donor count.
     pooled <- mbmm_pool_bank(
          particles,
          seed = 1
     )
     
     counts <- table(pooled$id)
     
     expect_true(all(counts == 5))
     expect_equal(nrow(pooled), 15)
     
     # A requested number below every donor's available count is retained.
     pooled_fixed <- mbmm_pool_bank(
          particles,
          n_per_donor = 3,
          seed = 1
     )
     
     expect_true(
          all(table(pooled_fixed$id) == 3)
     )
     
     expect_equal(
          nrow(pooled_fixed),
          9
     )
})


test_that("mbmm_pool_bank stops when equal donor contribution is impossible", {
     particles <- data.frame(
          id = rep(c("A", "B", "C"), times = c(10, 5, 8)),
          x = rnorm(23)
     )
     
     expect_error(
          mbmm_pool_bank(
               particles,
               n_per_donor = 6,
               seed = 1
          ),
          "fewer than 6 accepted particles"
     )
})


test_that("mbmm_pool_bank can intentionally allow fewer particles", {
     particles <- data.frame(
          id = rep(c("A", "B", "C"), times = c(10, 5, 8)),
          x = rnorm(23)
     )
     
     pooled <- mbmm_pool_bank(
          particles,
          n_per_donor = 6,
          allow_fewer = TRUE,
          seed = 1
     )
     
     expect_equal(
          as.integer(table(pooled$id)[c("A", "B", "C")]),
          c(6, 5, 6)
     )
})


test_that("mbmm_pool_bank applies eligible donor filtering", {
     particles <- data.frame(
          id = rep(c("A", "B", "C"), each = 5),
          x = rnorm(15)
     )
     
     pooled <- mbmm_pool_bank(
          particles,
          n_per_donor = 3,
          eligible_donors = c("A", "C"),
          seed = 1
     )
     
     expect_setequal(
          unique(pooled$id),
          c("A", "C")
     )
     
     expect_true(
          all(table(pooled$id) == 3)
     )
})


test_that("mbmm_join_targets builds complete donor-bundle-target combinations", {
     by_target <- list(
          IgG = data.frame(
               id = rep(c("A", "B"), each = 5),
               v = 1:10
          ),
          IgM = data.frame(
               id = rep(c("A", "B"), each = 5),
               v = 11:20
          )
     )
     
     bundles <- mbmm_join_targets(
          by_target,
          n_bundles = 4,
          seed = 1
     )
     
     expect_setequal(
          unique(bundles$target),
          c("IgG", "IgM")
     )
     
     # donors x bundles x targets
     expect_equal(
          nrow(bundles),
          2 * 4 * 2
     )
     
     # bundle_id runs from 1 to 4 within each donor.
     expect_setequal(
          unique(bundles$bundle_id),
          1:4
     )
     
     donor_bundle <- unique(
          bundles[c("id", "bundle_id")]
     )
     
     expect_equal(
          nrow(donor_bundle),
          2 * 4
     )
     
     # Every donor-bundle must contain both targets exactly once.
     bundle_target_counts <- table(
          bundles$id,
          bundles$bundle_id
     )
     
     expect_true(
          all(bundle_target_counts == 2)
     )
     
     expect_false(
          anyDuplicated(
               bundles[c("id", "bundle_id", "target")]
          ) > 0
     )
})


test_that("mbmm_join_targets uses only donors present in every target", {
     by_target <- list(
          IgG = data.frame(
               id = rep(c("A", "B"), each = 3),
               v = 1:6
          ),
          IgM = data.frame(
               id = rep(c("A", "C"), each = 3),
               v = 7:12
          )
     )
     
     bundles <- mbmm_join_targets(
          by_target,
          n_bundles = 2,
          seed = 1
     )
     
     expect_equal(
          unique(bundles$id),
          "A"
     )
     
     expect_equal(
          nrow(bundles),
          1 * 2 * 2
     )
})


test_that("mbmm_predict_bank evaluates simulate_fn at requested times", {
     bank <- data.frame(
          id = c("A", "B"),
          a = c(1, 2),
          b = c(0.5, -0.5)
     )
     
     pred <- mbmm_predict_bank(
          bank = bank,
          simulate_fn = linear_simulate,
          times = c(0, 1, 2)
     )
     
     expect_equal(
          nrow(pred),
          6
     )
     
     expect_equal(
          sort(unique(pred$time)),
          c(0, 1, 2)
     )
     
     row_a0 <- pred[
          pred$id == "A" &
               pred$time == 0,
          ,
          drop = FALSE
     ]
     
     expect_equal(
          row_a0$pred,
          1
     )
     
     expect_equal(
          attr(pred, "n_failed_bank_rows"),
          0
     )
})


test_that("mbmm_predict_bank detects invalid simulation output", {
     bank <- data.frame(
          id = "A",
          a = 1,
          b = 0.5
     )
     
     bad_simulate <- function(theta_row, times) {
          c(1, 2)
     }
     
     expect_error(
          mbmm_predict_bank(
               bank = bank,
               simulate_fn = bad_simulate,
               times = c(0, 1, 2),
               failure = "error"
          ),
          "Expected 3 finite predictions"
     )
})

test_that("mbmm_abc_fit supports acceptance-fraction retention", {
     times <- 0:3
     
     obs <- data.frame(
          time = times,
          y = 1.5 + 0.8 * times
     )
     
     fit <- mbmm_abc_fit(
          id = "S1",
          obs = obs,
          simulate_fn = linear_simulate,
          prior_sample_fn = linear_prior,
          n_draws = 1000,
          n_keep = NULL,
          accept_frac = 0.025,
          seed = 1
     )
     
     expect_equal(nrow(fit), 25)
     expect_equal(attr(fit, "n_keep"), 25)
     expect_equal(attr(fit, "n_draws"), 1000)
     expect_equal(
          attr(fit, "acceptance_fraction_total"),
          0.025
     )
})


test_that("mbmm_abc_fit discards failed simulations and reports QC", {
     simulate_with_failures <- function(theta_row, times) {
          if (theta_row$a < 2.5) {
               stop("intentional test failure")
          }
          
          as.numeric(theta_row$a) +
               as.numeric(theta_row$b) * times
     }
     
     obs <- data.frame(
          time = 0:2,
          y = c(3, 3.5, 4)
     )
     
     fit <- mbmm_abc_fit(
          id = "S1",
          obs = obs,
          simulate_fn = simulate_with_failures,
          prior_sample_fn = linear_prior,
          n_draws = 500,
          n_keep = 20,
          seed = 2,
          failure = "discard"
     )
     
     expect_equal(nrow(fit), 20)
     expect_gt(attr(fit, "n_failed"), 0)
     expect_equal(
          attr(fit, "n_valid") +
               attr(fit, "n_failed"),
          500
     )
})