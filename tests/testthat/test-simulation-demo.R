test_that("synthetic antibody generator returns valid standardized values", {
     x <- mbmm_simulate_demo_antibody(
          n = 100,
          seed = 101
     )
     
     expect_s3_class(x, "data.frame")
     expect_equal(nrow(x), 100L)
     
     expect_true(
          all(
               c(
                    "antibody",
                    "log_antibody",
                    "z_antibody",
                    "high_component",
                    "calibration_mean",
                    "calibration_sd"
               ) %in% names(x)
          )
     )
     
     expect_true(all(is.finite(x$antibody)))
     expect_true(all(x$antibody > 0))
     expect_true(all(is.finite(x$log_antibody)))
     expect_true(all(is.finite(x$z_antibody)))
     
     expect_equal(mean(x$z_antibody), 0, tolerance = 1e-12)
     expect_equal(stats::sd(x$z_antibody), 1, tolerance = 1e-12)
     
     expect_true(all(x$high_component %in% c(0L, 1L)))
})


test_that("synthetic antibody generation is reproducible", {
     x1 <- mbmm_simulate_demo_antibody(
          n = 50,
          seed = 202
     )
     
     x2 <- mbmm_simulate_demo_antibody(
          n = 50,
          seed = 202
     )
     
     expect_identical(x1, x2)
})


test_that("public visit schedules are valid", {
     schedules <- mbmm_demo_visit_schedules()
     
     expect_type(schedules, "list")
     expect_true("complete" %in% names(schedules))
     expect_true(length(schedules) >= 2L)
     
     valid <- vapply(
          schedules,
          function(x) {
               is.numeric(x) &&
                    length(x) >= 2L &&
                    all(is.finite(x)) &&
                    all(x >= 0) &&
                    identical(x, sort(unique(x)))
          },
          logical(1)
     )
     
     expect_true(all(valid))
})


test_that("each transient or abortive subtype can be generated", {
     times <- c(0, 2, 4, 7, 10, 14, 21)
     subtypes <- c(
          "noise_like",
          "single_blip",
          "two_point",
          "fast_rise_drop"
     )
     
     results <- lapply(
          seq_along(subtypes),
          function(i) {
               mbmm_simulate_abortive_pattern(
                    times = times,
                    lod = 100,
                    sigma_obs = 0.3,
                    subtype = subtypes[i],
                    seed = 300 + i
               )
          }
     )
     
     for (i in seq_along(results)) {
          x <- results[[i]]
          
          expect_s3_class(x, "data.frame")
          expect_equal(nrow(x), length(times))
          expect_equal(unique(x$subtype), subtypes[i])
          
          expect_true(
               all(
                    c(
                         "time",
                         "y",
                         "censored",
                         "lod",
                         "log10_true",
                         "log10_observed",
                         "subtype"
                    ) %in% names(x)
               )
          )
          
          expect_true(all(is.finite(x$time)))
          expect_true(all(is.finite(x$y)))
          expect_true(all(x$y > 0))
          expect_true(all(x$lod == 100))
          expect_true(all(x$censored %in% c(TRUE, FALSE)))
     }
})


test_that("intended abortive patterns have the expected latent shape", {
     times <- c(0, 2, 4, 7, 10, 14, 21)
     lod <- 100
     lod_log10 <- log10(lod)
     
     single <- mbmm_simulate_abortive_pattern(
          times = times,
          lod = lod,
          sigma_obs = 0,
          subtype = "single_blip",
          seed = 401
     )
     
     two <- mbmm_simulate_abortive_pattern(
          times = times,
          lod = lod,
          sigma_obs = 0,
          subtype = "two_point",
          seed = 402
     )
     
     fast <- mbmm_simulate_abortive_pattern(
          times = times,
          lod = lod,
          sigma_obs = 0,
          subtype = "fast_rise_drop",
          seed = 403
     )
     
     expect_equal(sum(single$log10_true >= lod_log10), 1L)
     expect_equal(sum(two$log10_true >= lod_log10), 2L)
     expect_equal(sum(fast$log10_true >= lod_log10), 1L)
})


test_that("manuscript-style simulation returns aligned package inputs", {
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 40,
          prevalence = 0.30,
          beta_true = -1.5,
          p_abortive = 0.40,
          n_bank_donors = 6,
          seed = 501
     )
     
     expect_type(demo, "list")
     
     expect_true(
          all(
               c(
                    "subject_data",
                    "truth",
                    "observations",
                    "obs_list",
                    "covariates",
                    "productive_parameter_bank",
                    "settings"
               ) %in% names(demo)
          )
     )
     
     expect_equal(nrow(demo$subject_data), 40L)
     expect_equal(nrow(demo$truth), 40L)
     expect_equal(length(demo$obs_list), 40L)
     expect_equal(nrow(demo$covariates), 40L)
     
     expect_identical(
          names(demo$obs_list),
          demo$subject_data$id
     )
     
     expect_identical(
          demo$covariates$id,
          demo$subject_data$id
     )
     
     expect_equal(
          demo$covariates$baseline_antibody,
          demo$subject_data$baseline_antibody
     )
     
     required_obs_columns <- c(
          "time",
          "y",
          "censored",
          "lod"
     )
     
     expect_true(
          all(
               vapply(
                    demo$obs_list,
                    function(x) {
                         is.data.frame(x) &&
                              all(required_obs_columns %in% names(x)) &&
                              all(is.finite(x$time)) &&
                              all(is.finite(x$y)) &&
                              all(x$y > 0) &&
                              all(is.finite(x$lod)) &&
                              all(x$lod > 0)
                    },
                    logical(1)
               )
          )
     )
})


test_that("marginal productive probability is calibrated to the target", {
     target <- 0.30
     
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 60,
          prevalence = target,
          beta_true = -1.5,
          p_abortive = 0.20,
          n_bank_donors = 6,
          seed = 601
     )
     
     expect_equal(
          mean(demo$subject_data$productive_probability),
          target,
          tolerance = 1e-4
     )
     
     expect_equal(
          demo$settings$realized_mean_probability,
          target,
          tolerance = 1e-4
     )
})


test_that("negative biomarker effect lowers productive-infection probability", {
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 80,
          prevalence = 0.30,
          beta_true = -1.5,
          p_abortive = 0.20,
          n_bank_donors = 6,
          seed = 701
     )
     
     association <- stats::cor(
          demo$subject_data$baseline_antibody,
          demo$subject_data$productive_probability
     )
     
     expect_true(is.finite(association))
     expect_lt(association, 0)
})


test_that("bank-eligible subjects are productive and receive the complete schedule", {
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 50,
          prevalence = 0.35,
          beta_true = -1.2,
          p_abortive = 0.40,
          n_bank_donors = 8,
          seed = 801
     )
     
     eligible <- demo$subject_data$bank_eligible
     
     expect_true(any(eligible))
     expect_true(all(demo$subject_data$productive_true[eligible] == 1L))
     expect_true(
          all(
               demo$subject_data$visit_schedule[eligible] == "complete"
          )
     )
     
     complete_times <- mbmm_demo_visit_schedules()[["complete"]]
     
     expect_true(
          all(
               vapply(
                    demo$obs_list[demo$subject_data$id[eligible]],
                    function(x) identical(x$time, complete_times),
                    logical(1)
               )
          )
     )
})


test_that("productive trajectories satisfy the requested minimum detection rule", {
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 40,
          prevalence = 0.35,
          beta_true = -1.0,
          p_abortive = 0.20,
          n_bank_donors = 5,
          productive_min_detect = 1,
          productive_max_tries = 50,
          seed = 901
     )
     
     productive <- demo$subject_data$productive_true == 1L
     
     expect_true(any(productive))
     expect_true(all(demo$subject_data$n_detect[productive] >= 1L))
})


test_that("nonestablishment subjects have defined abortive labels", {
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 40,
          prevalence = 0.30,
          beta_true = -1.0,
          p_abortive = 0.40,
          n_bank_donors = 5,
          seed = 1001
     )
     
     nonestablishment <- demo$subject_data$productive_true == 0L
     productive <- demo$subject_data$productive_true == 1L
     
     allowed <- c(
          "noise_like",
          "single_blip",
          "two_point",
          "fast_rise_drop"
     )
     
     expect_true(
          all(
               demo$subject_data$abortive_subtype[nonestablishment] %in% allowed
          )
     )
     
     expect_true(
          all(
               is.na(demo$subject_data$abortive_subtype[productive])
          )
     )
})


test_that("complete dataset simulation is reproducible", {
     demo1 <- mbmm_simulate_manuscript_demo(
          n_subjects = 30,
          prevalence = 0.30,
          beta_true = -1.2,
          p_abortive = 0.40,
          n_bank_donors = 4,
          seed = 1101
     )
     
     demo2 <- mbmm_simulate_manuscript_demo(
          n_subjects = 30,
          prevalence = 0.30,
          beta_true = -1.2,
          p_abortive = 0.40,
          n_bank_donors = 4,
          seed = 1101
     )
     
     expect_identical(demo1$subject_data, demo2$subject_data)
     expect_identical(demo1$observations, demo2$observations)
     expect_identical(demo1$obs_list, demo2$obs_list)
     expect_identical(demo1$settings, demo2$settings)
})


test_that("reduced scenario grid has the expected structure", {
     grid <- mbmm_demo_scenario_grid(
          replicates = 2,
          n_subjects = 50,
          prevalence = 0.30,
          beta_grid = c(0, -1.5),
          p_abortive_grid = c(0, 0.40),
          base_seed = 20260724
     )
     
     expect_s3_class(grid, "data.frame")
     expect_equal(nrow(grid), 8L)
     
     expect_true(
          all(
               c(
                    "n_subjects",
                    "prevalence",
                    "beta_true",
                    "p_abortive",
                    "replicate",
                    "scenario",
                    "seed"
               ) %in% names(grid)
          )
     )
     
     expect_equal(sort(unique(grid$beta_true)), c(-1.5, 0))
     expect_equal(sort(unique(grid$p_abortive)), c(0, 0.40))
     expect_equal(sort(unique(grid$replicate)), 1:2)
     expect_equal(length(unique(grid$seed)), nrow(grid))
})


test_that("simulation helpers reject invalid inputs", {
     expect_error(
          mbmm_simulate_demo_antibody(n = 1),
          "at least 2"
     )
     
     expect_error(
          mbmm_simulate_abortive_pattern(
               times = numeric(0)
          ),
          "`times`"
     )
     
     expect_error(
          mbmm_simulate_abortive_pattern(
               times = c(0, 2),
               lod = 0
          ),
          "`lod`"
     )
     
     expect_error(
          mbmm_simulate_manuscript_demo(
               n_subjects = 5
          ),
          "at least 10"
     )
     
     expect_error(
          mbmm_simulate_manuscript_demo(
               prevalence = 1
          ),
          "`prevalence`"
     )
     
     expect_error(
          mbmm_demo_scenario_grid(
               replicates = 0
          ),
          "`replicates`"
     )
})