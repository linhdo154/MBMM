#' @noRd
.mbmm_run_demo_scenario <- function(
          setting,
          productive_bank,
          sigma = 0.50,
          lod_fraction = 0.5,
          ridge = 1e-4
) {
     if (!is.data.frame(setting) || nrow(setting) != 1L) {
          stop("`setting` must be a one-row data frame.", call. = FALSE)
     }
     
     demo <- mbmm_simulate_manuscript_demo(
          n_subjects = setting$n_subjects,
          prevalence = setting$prevalence,
          beta_true = setting$beta_true,
          p_abortive = setting$p_abortive,
          n_bank_donors = 0L,
          sigma_obs = sigma,
          lod = 100,
          seed = setting$seed
     )
     
     fit <- mbmm_productive_infection(
          obs_list = demo$obs_list,
          productive_bank = productive_bank,
          sigma = sigma,
          formula = ~baseline_antibody,
          data = demo$covariates,
          productive_class = "productive",
          nonestablishment_class = "nonestablishment",
          lod_fraction = lod_fraction,
          penalty = "ridge",
          ridge = ridge,
          id_col = "id",
          donor_col = "id",
          bundle_col = "bundle_id",
          log_scale = FALSE,
          donor_weighting = "equal",
          loo_productive = FALSE,
          verbose = FALSE
     )
     
     list(
          setting = setting,
          demo = demo,
          fit = fit
     )
}

#' Run a manuscript-style simulation grid
#'
#' Runs independent simulation scenarios sequentially or across multiple
#' worker processes.
#'
#' @param grid A simulation grid returned by [mbmm_demo_scenario_grid()].
#' @param productive_bank Predicted productive reference bank.
#' @param n_cores Number of worker processes. Use `1` for sequential execution.
#'
#' @return A list containing one fitted result for each row of `grid`.
#' @export
mbmm_run_simulation_grid <- function(
          grid,
          productive_bank,
          n_cores = 1L
) {
     if (!is.data.frame(grid) || nrow(grid) < 1L) {
          stop("`grid` must be a nonempty data frame.", call. = FALSE)
     }
     
     if (
          length(n_cores) != 1L ||
          !is.numeric(n_cores) ||
          !is.finite(n_cores) ||
          n_cores < 1
     ) {
          stop("`n_cores` must be at least 1.", call. = FALSE)
     }
     
     n_cores <- min(
          as.integer(n_cores),
          nrow(grid)
     )
     
     grid_rows <- lapply(
          seq_len(nrow(grid)),
          function(i) {
               grid[i, , drop = FALSE]
          }
     )
     
     if (n_cores == 1L) {
          return(
               lapply(
                    grid_rows,
                    .mbmm_run_demo_scenario,
                    productive_bank = productive_bank
               )
          )
     }
     
     cluster <- parallel::makeCluster(n_cores)
     
     on.exit(
          parallel::stopCluster(cluster),
          add = TRUE
     )
     
     parallel::clusterEvalQ(
          cluster,
          library(mbmm)
     )
     
     parallel::parLapply(
          cluster,
          grid_rows,
          function(setting, productive_bank) {
               worker <- utils::getFromNamespace(
                    ".mbmm_run_demo_scenario",
                    "mbmm"
               )
               
               worker(
                    setting = setting,
                    productive_bank = productive_bank
               )
          },
          productive_bank = productive_bank
     )
}