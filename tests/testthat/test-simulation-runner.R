test_that("simulation grid runner works sequentially", {
     reference_demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 35,
          prevalence = 0.40,
          beta_true = 0,
          p_abortive = 0,
          n_bank_donors = 4,
          sigma_obs = 0.40,
          lod = 100,
          seed = 20260730
     )
     
     reference_ids <- reference_demo$subject_data$id[
          reference_demo$subject_data$bank_eligible
     ]
     
     abc_particles <- mbmm_abc_bank(
          obs_list = reference_demo$obs_list[reference_ids],
          simulate_fn = mbmm_example_eclipse_simulate,
          prior_sample_fn = mbmm_example_eclipse_prior,
          n_draws = 100,
          n_keep = 10,
          seed = 20260731,
          verbose = FALSE
     )
     
     pooled_productive <- mbmm_pool_bank(
          particles = abc_particles,
          seed = 20260801
     )
     
     prediction_times <- sort(
          unique(
               unlist(
                    mbmm_demo_visit_schedules(),
                    use.names = FALSE
               )
          )
     )
     
     productive_bank <- mbmm_predict_bank(
          bank = pooled_productive,
          simulate_fn = mbmm_example_eclipse_simulate,
          times = prediction_times
     )
     
     grid <- mbmm_demo_scenario_grid(
          replicates = 2,
          n_subjects = 30,
          prevalence = 0.30,
          beta_grid = -1.5,
          p_abortive_grid = 0.40,
          base_seed = 123
     )
     
     runs <- mbmm_run_simulation_grid(
          grid = grid,
          productive_bank = productive_bank,
          n_cores = 1
     )
     
     expect_type(runs, "list")
     expect_length(runs, 2)
     
     expect_true(
          all(
               vapply(
                    runs,
                    function(x) {
                         identical(
                              names(x),
                              c("setting", "demo", "fit")
                         )
                    },
                    logical(1)
               )
          )
     )
     
     expect_equal(
          vapply(
               runs,
               function(x) x$setting$seed,
               numeric(1)
          ),
          c(123, 124)
     )
})


test_that("parallel and sequential simulation runs agree", {
     skip_on_cran()
     
     n_available <- parallel::detectCores(
          logical = FALSE
     )
     
     skip_if(
          is.na(n_available) || n_available < 2,
          "At least two physical cores are required."
     )
     
     reference_demo <- mbmm_simulate_manuscript_demo(
          n_subjects = 35,
          prevalence = 0.40,
          beta_true = 0,
          p_abortive = 0,
          n_bank_donors = 4,
          sigma_obs = 0.40,
          lod = 100,
          seed = 20260802
     )
     
     reference_ids <- reference_demo$subject_data$id[
          reference_demo$subject_data$bank_eligible
     ]
     
     abc_particles <- mbmm_abc_bank(
          obs_list = reference_demo$obs_list[reference_ids],
          simulate_fn = mbmm_example_eclipse_simulate,
          prior_sample_fn = mbmm_example_eclipse_prior,
          n_draws = 100,
          n_keep = 10,
          seed = 20260803,
          verbose = FALSE
     )
     
     pooled_productive <- mbmm_pool_bank(
          particles = abc_particles,
          seed = 20260804
     )
     
     prediction_times <- sort(
          unique(
               unlist(
                    mbmm_demo_visit_schedules(),
                    use.names = FALSE
               )
          )
     )
     
     productive_bank <- mbmm_predict_bank(
          bank = pooled_productive,
          simulate_fn = mbmm_example_eclipse_simulate,
          times = prediction_times
     )
     
     grid <- mbmm_demo_scenario_grid(
          replicates = 2,
          n_subjects = 30,
          prevalence = 0.30,
          beta_grid = -1.5,
          p_abortive_grid = 0.40,
          base_seed = 456
     )
     
     sequential <- mbmm_run_simulation_grid(
          grid = grid,
          productive_bank = productive_bank,
          n_cores = 1
     )
     
     parallel_runs <- mbmm_run_simulation_grid(
          grid = grid,
          productive_bank = productive_bank,
          n_cores = 2
     )
     
     sequential_seeds <- vapply(
          sequential,
          function(x) x$setting$seed,
          numeric(1)
     )
     
     parallel_seeds <- vapply(
          parallel_runs,
          function(x) x$setting$seed,
          numeric(1)
     )
     
     sequential_lrt <- vapply(
          sequential,
          function(x) x$fit$gate$LRT,
          numeric(1)
     )
     
     parallel_lrt <- vapply(
          parallel_runs,
          function(x) x$fit$gate$LRT,
          numeric(1)
     )
     
     expect_identical(
          sequential_seeds,
          parallel_seeds
     )
     
     expect_equal(
          sequential_lrt,
          parallel_lrt,
          tolerance = 1e-10
     )
})