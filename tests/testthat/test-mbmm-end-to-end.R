# End-to-end pipeline test using a fast closed-form trajectory model
# (single-exponential decay) instead of the deSolve-based eclipse example,
# so this runs quickly on CRAN.
decay_simulate <- function(theta_row, times) {
     as.numeric(theta_row$A0) *
          exp(-as.numeric(theta_row$k) * times)
}

decay_prior <- function(n) {
     data.frame(
          A0 = 10^stats::runif(n, 0.5, 1.5),
          k = stats::runif(n, 0.05, 1.0)
     )
}


test_that("the full mbmm() pipeline runs end to end and separates two truly different groups", {
     set.seed(11)
     
     times <- c(0, 1, 2, 4, 7)
     
     # Fast-decay truth versus slow-decay truth.
     make_subject <- function(id, k_true, x) {
          y <- 10 *
               exp(-k_true * times) *
               exp(rnorm(length(times), 0, 0.05))
          
          list(
               obs = data.frame(
                    time = times,
                    y = y
               ),
               x = x
          )
     }
     
     n_per_group <- 12
     
     fast_x <- rnorm(
          n_per_group,
          mean = 1
     )
     
     slow_x <- rnorm(
          n_per_group,
          mean = -1
     )
     
     subs <- c(
          Map(
               function(i, x) {
                    make_subject(
                         paste0("F", i),
                         0.8,
                         x
                    )
               },
               seq_len(n_per_group),
               fast_x
          ),
          Map(
               function(i, x) {
                    make_subject(
                         paste0("S", i),
                         0.1,
                         x
                    )
               },
               seq_len(n_per_group),
               slow_x
          )
     )
     
     ids <- c(
          paste0("F", seq_len(n_per_group)),
          paste0("S", seq_len(n_per_group))
     )
     
     names(subs) <- NULL
     
     obs_list <- stats::setNames(
          lapply(subs, `[[`, "obs"),
          ids
     )
     
     covariates <- data.frame(
          id = ids,
          xcov = vapply(
               subs,
               `[[`,
               numeric(1),
               "x"
          )
     )
     
     # Build ABC banks from a subset of donors in each group.
     fast_donors <- paste0("F", 1:6)
     slow_donors <- paste0("S", 1:6)
     
     bank_fast <- mbmm_abc_bank(
          obs_list = obs_list[fast_donors],
          simulate_fn = decay_simulate,
          prior_sample_fn = decay_prior,
          n_draws = 3000,
          n_keep = 40,
          seed = 1,
          verbose = FALSE
     )
     
     bank_slow <- mbmm_abc_bank(
          obs_list = obs_list[slow_donors],
          simulate_fn = decay_simulate,
          prior_sample_fn = decay_prior,
          n_draws = 3000,
          n_keep = 40,
          seed = 2,
          verbose = FALSE
     )
     
     pooled_fast <- mbmm_pool_bank(
          particles = bank_fast,
          seed = 1
     )
     
     pooled_slow <- mbmm_pool_bank(
          particles = bank_slow,
          seed = 1
     )
     
     pred_fast <- mbmm_predict_bank(
          bank = pooled_fast,
          simulate_fn = decay_simulate,
          times = times
     )
     
     pred_slow <- mbmm_predict_bank(
          bank = pooled_slow,
          simulate_fn = decay_simulate,
          times = times
     )
     
     fit <- mbmm(
          obs_list = obs_list,
          bank_pred_list = list(
               fast = pred_fast,
               slow = pred_slow
          ),
          sigma = 0.15,
          formula = ~xcov,
          data = covariates,
          classes = c("fast", "slow"),
          id_col = "id",
          log_scale = FALSE,
          donor_weighting = "equal",
          loo = TRUE,
          penalty = "ridge",
          ridge = 1e-4,
          verbose = FALSE
     )
     
     expect_s3_class(
          fit,
          "mbmm_fit"
     )
     
     expect_equal(
          nrow(fit$logL),
          length(ids)
     )
     
     expect_equal(
          fit$classes,
          c("fast", "slow")
     )
     
     expect_equal(
          fit$scoring$type,
          "bank-versus-bank"
     )
     
     expect_true(
          isTRUE(fit$scoring$loo)
     )
     
     expect_equal(
          fit$scoring$donor_weighting,
          "equal"
     )
     
     expect_false(
          fit$scoring$log_scale
     )
     
     expect_equal(
          fit$gate$convergence_full,
          0
     )
     
     expect_equal(
          fit$gate$convergence_null,
          0
     )
     
     # Higher xcov should indicate a higher probability of fast decay.
     beta_hat <- fit$gate$coef$estimate[
          fit$gate$coef$term == "xcov"
     ]
     
     expect_gt(
          beta_hat,
          0
     )
     
     # MAP class recovery should be well above chance.
     truth <- ifelse(
          startsWith(fit$logL$id, "F"),
          "fast",
          "slow"
     )
     
     accuracy <- mean(
          fit$gate$posterior$map_class == truth
     )
     
     expect_gt(
          accuracy,
          0.75
     )
     
     expect_output(
          print(fit),
          "MBMM fit"
     )
     
     summary_result <- summary(fit)
     
     expect_s3_class(
          summary_result,
          "mbmm_fit"
     )
})