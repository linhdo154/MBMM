test_that("mbmm_backtransform_threshold handles common biomarker scales", {
     expect_equal(
          mbmm_backtransform_threshold(
               threshold = 1,
               center = 2,
               scale = 0.5,
               transform = "identity"
          ),
          2.5
     )
     
     expect_equal(
          mbmm_backtransform_threshold(
               threshold = 1,
               center = 2,
               scale = 0.5,
               transform = "log10"
          ),
          10^2.5
     )
     
     expect_equal(
          mbmm_backtransform_threshold(
               threshold = 1,
               center = 2,
               scale = 0.5,
               transform = "log"
          ),
          exp(2.5)
     )
})


test_that("bootstrap gate threshold returns estimates and diagnostics", {
     set.seed(42)
     
     n <- 120
     x <- stats::rnorm(n)
     
     z <- stats::rbinom(
          n,
          size = 1,
          prob = mbmm_inv_logit(-0.5 + 1.5 * x)
     )
     
     logL <- data.frame(
          id = paste0("S", seq_len(n)),
          logL_pos = ifelse(
               z == 1,
               log(0.98),
               log(0.02)
          ),
          logL_neg = ifelse(
               z == 1,
               log(0.02),
               log(0.98)
          )
     )
     
     data <- data.frame(
          id = logL$id,
          x = x
     )
     
     result <- mbmm_bootstrap_gate_threshold(
          logL = logL,
          formula = ~x,
          data = data,
          focal = "x",
          classes = c("pos", "neg"),
          n_boot = 30,
          seed = 1,
          verbose = FALSE
     )
     
     expect_s3_class(
          result,
          "mbmm_threshold_summary"
     )
     
     expect_true(
          is.finite(result$point$threshold_model)
     )
     
     expect_length(
          result$confidence_interval$model_scale,
          2
     )
     
     expect_length(
          result$confidence_interval$original_scale,
          2
     )
     
     expect_equal(
          nrow(result$bootstrap),
          30
     )
     
     expect_gt(
          result$diagnostics$n_threshold_success,
          0
     )
     
     expect_true(
          result$pct_extrapolated >= 0 &&
               result$pct_extrapolated <= 100
     )
     
     expect_output(
          print(result),
          "MBMM gate-threshold summary"
     )
     
     expect_s3_class(
          summary(result),
          "mbmm_threshold_summary"
     )
})


test_that("bootstrap gate threshold supports log10 back-transformation", {
     set.seed(43)
     
     n <- 100
     raw_biomarker <- 10^stats::rnorm(
          n,
          mean = 2,
          sd = 0.4
     )
     
     log_biomarker <- log10(
          raw_biomarker
     )
     
     center <- mean(
          log_biomarker
     )
     
     scale <- stats::sd(
          log_biomarker
     )
     
     biomarker_z <- (
          log_biomarker -
               center
     ) / scale
     
     z <- stats::rbinom(
          n,
          size = 1,
          prob = mbmm_inv_logit(
               0.3 -
                    1.2 * biomarker_z
          )
     )
     
     logL <- data.frame(
          id = paste0("S", seq_len(n)),
          logL_productive = ifelse(
               z == 1,
               log(0.98),
               log(0.02)
          ),
          logL_nonestablishment = ifelse(
               z == 1,
               log(0.02),
               log(0.98)
          )
     )
     
     data <- data.frame(
          id = logL$id,
          biomarker_z = biomarker_z
     )
     
     result <- mbmm_bootstrap_gate_threshold(
          logL = logL,
          formula = ~biomarker_z,
          data = data,
          focal = "biomarker_z",
          classes = c(
               "productive",
               "nonestablishment"
          ),
          center = center,
          scale = scale,
          transform = "log10",
          n_boot = 20,
          seed = 2,
          verbose = FALSE
     )
     
     expected_original <- 10^(
          center +
               scale *
               result$point$threshold_model
     )
     
     expect_equal(
          result$point$threshold_original,
          expected_original
     )
})