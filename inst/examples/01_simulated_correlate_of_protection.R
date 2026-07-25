# Standalone script version of vignette("mbmm-intro"). Run with:
#   source(system.file("examples", "01_simulated_correlate_of_protection.R", package = "mbmm"))
library(mbmm)

set.seed(20260722)

## 1. Simulate a demo correlate-of-protection cohort ------------------------
demo <- mbmm_simulate_cop_demo(
  n_subjects = 70, prevalence = 0.4, beta_true = -1.3,
  n_bank_donors = 10, sigma_obs = 0.3, lod = 3, seed = 20260722
)

cat("n subjects:", length(demo$obs_list), "\n")
cat("n truly infected (Z=1):", sum(demo$truth$Z_true), "\n")
cat("n bank-eligible donors:", sum(demo$truth$bank_eligible), "\n")

bank_ids <- demo$truth$id[demo$truth$bank_eligible]
str(demo$obs_list[[bank_ids[1]]])

## 2. ABC-fit the "productive infection" reference bank ---------------------
t0 <- Sys.time()
abc_particles <- mbmm_abc_bank(
  obs_list = demo$obs_list[bank_ids],
  simulate_fn = mbmm_example_eclipse_simulate,
  prior_sample_fn = mbmm_example_eclipse_prior,
  n_draws = 600, n_keep = 60, seed = 1, verbose = FALSE
)
cat("ABC time:", Sys.time() - t0, "\n")
print(head(abc_particles))

pooled_infected <- mbmm_pool_bank(abc_particles, seed = 1)
cat("pooled bank rows:", nrow(pooled_infected), "\n")

visit_grid <- c(0, 2, 4, 6, 9, 12, 16, 21)
pred_infected <- mbmm_predict_bank(pooled_infected, mbmm_example_eclipse_simulate,
                                    times = visit_grid)
cat("pred_infected rows:", nrow(pred_infected), "\n")

## 3. Build the "non-infected" reference bank (flat near-LOD trajectory) ----
pred_noninfected <- data.frame(
  id = "reference", bundle_id = 1,
  time = visit_grid, pred = rep(1.5, length(visit_grid))
)

## 4. Score every subject against both banks (log scale) --------------------
obs_log <- lapply(demo$obs_list, function(d) {
  data.frame(time = d$time, y = mbmm_safe_log10(d$y),
             censored = d$censored, lod = mbmm_safe_log10(d$lod))
})
pred_infected_log <- pred_infected
pred_infected_log$pred <- mbmm_safe_log10(pred_infected$pred)
pred_noninfected_log <- pred_noninfected
pred_noninfected_log$pred <- mbmm_safe_log10(pred_noninfected$pred)

fit <- mbmm(
  obs_list = obs_log,
  bank_pred_list = list(infected = pred_infected_log, noninfected = pred_noninfected_log),
  sigma = 0.3,
  formula = ~baseline_antibody,
  data = demo$covariates,
  classes = c("infected", "noninfected"),
  ridge = 1e-4, n_perm = 200, seed = 1, id_col = "id", verbose = FALSE
)

print(fit)

## 5. Recovery check against ground truth -----------------------------------
truth <- demo$truth$Z_true[match(fit$logL$id, demo$truth$id)]
pred_map <- ifelse(fit$gate$posterior$map_class == "infected", 1, 0)
cat("\nConfusion matrix (rows = truth, cols = MBMM MAP class):\n")
print(table(truth, pred_map))
cat("Accuracy:", mean(truth == pred_map), "\n")

cat("\nTrue gate: alpha =", demo$gate_truth$alpha, " beta =", demo$gate_truth$beta, "\n")
cat("Fitted gate coefficients:\n")
print(fit$gate$coef)

thr <- mbmm_gate_threshold(fit$gate, "baseline_antibody")
cat("\nFitted 50% gate threshold (baseline_antibody):", thr, "\n")
