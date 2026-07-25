# Package index

## Fit the model

Top-level fitting functions

- [`mbmm()`](https://linhdo154.github.io/MBMM/reference/mbmm.md) : Fit a
  bank-based MBMM analysis
- [`mbmm_fit_gate()`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md)
  : Fit the MBMM gate model
- [`mbmm_abc_fit()`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md)
  : Fit an ABC rejection-sampling posterior for one subject

## ABC reference banks

Build and pool Approximate Bayesian Computation banks

- [`mbmm_abc_bank()`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_bank.md)
  : Fit ABC posteriors for multiple subjects or donors
- [`mbmm_pool_bank()`](https://linhdo154.github.io/MBMM/reference/mbmm_pool_bank.md)
  : Pool per-donor ABC particles into a mechanistic reference bank
- [`mbmm_predict_bank()`](https://linhdo154.github.io/MBMM/reference/mbmm_predict_bank.md)
  : Predict mechanistic bank trajectories at required observation times
- [`mbmm_join_targets()`](https://linhdo154.github.io/MBMM/reference/mbmm_join_targets.md)
  : Join target-specific particle banks into complete multitarget
  bundles

## Gate model

Covariate-linked latent-class gate

- [`mbmm_gate_threshold()`](https://linhdo154.github.io/MBMM/reference/mbmm_gate_threshold.md)
  : Fitted MBMM gate threshold for a continuous covariate
- [`mbmm_bootstrap_gate_threshold()`](https://linhdo154.github.io/MBMM/reference/mbmm_bootstrap_gate_threshold.md)
  : Bootstrap an MBMM gate threshold
- [`mbmm_backtransform_threshold()`](https://linhdo154.github.io/MBMM/reference/mbmm_backtransform_threshold.md)
  : Back-transform an MBMM gate threshold
- [`mbmm_predict_gate()`](https://linhdo154.github.io/MBMM/reference/mbmm_predict_gate.md)
  : Predict fitted MBMM gate probabilities
- [`summary(`*`<mbmm_gate>`*`)`](https://linhdo154.github.io/MBMM/reference/summary.mbmm_gate.md)
  : Summarize an mbmm_gate object
- [`print(`*`<mbmm_gate>`*`)`](https://linhdo154.github.io/MBMM/reference/print.mbmm_gate.md)
  : Print an mbmm_gate object
- [`plot(`*`<mbmm_gate>`*`)`](https://linhdo154.github.io/MBMM/reference/plot.mbmm_gate.md)
  : Plot a fitted MBMM gate curve against one continuous covariate

## Scoring and thresholds

- [`mbmm_score_subjects()`](https://linhdo154.github.io/MBMM/reference/mbmm_score_subjects.md)
  : Score every subject against class-specific reference banks
- [`mbmm_class_loglik()`](https://linhdo154.github.io/MBMM/reference/mbmm_class_loglik.md)
  : Reference-bank class-conditional log-likelihood for one subject
- [`mbmm_noise_floor_loglik()`](https://linhdo154.github.io/MBMM/reference/mbmm_noise_floor_loglik.md)
  : Analytic non-establishment class likelihood for COVID5001-style
  analyses
- [`summary(`*`<mbmm_threshold_summary>`*`)`](https://linhdo154.github.io/MBMM/reference/summary.mbmm_threshold_summary.md)
  : Summarize an MBMM threshold summary
- [`print(`*`<mbmm_threshold_summary>`*`)`](https://linhdo154.github.io/MBMM/reference/print.mbmm_threshold_summary.md)
  : Print an MBMM threshold summary

## Simulation and demo data

Simulate trajectories and example datasets used in the vignettes

- [`mbmm_simulate_demo_antibody()`](https://linhdo154.github.io/MBMM/reference/mbmm_simulate_demo_antibody.md)
  : Simulate a fully synthetic baseline antibody measurement
- [`mbmm_simulate_abortive_pattern()`](https://linhdo154.github.io/MBMM/reference/mbmm_simulate_abortive_pattern.md)
  : Simulate one transient or abortive detection pattern
- [`mbmm_simulate_cop_demo()`](https://linhdo154.github.io/MBMM/reference/mbmm_simulate_cop_demo.md)
  : Simulate a synthetic correlate-of-protection cohort
- [`mbmm_simulate_manuscript_demo()`](https://linhdo154.github.io/MBMM/reference/mbmm_simulate_manuscript_demo.md)
  : Simulate a reduced manuscript-style MBMM dataset
- [`mbmm_example_eclipse_simulate()`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_simulate.md)
  : Example mechanistic trajectory model: eclipse-phase viral kinetics
- [`mbmm_example_eclipse_prior()`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_prior.md)
  : Example ABC prior sampler for the eclipse-phase model
- [`mbmm_demo_visit_schedules()`](https://linhdo154.github.io/MBMM/reference/mbmm_demo_visit_schedules.md)
  : Public visit schedules for the reduced MBMM simulation
- [`mbmm_demo_scenario_grid()`](https://linhdo154.github.io/MBMM/reference/mbmm_demo_scenario_grid.md)
  : Create the default reduced manuscript-style scenario grid
- [`mbmm_run_simulation_grid()`](https://linhdo154.github.io/MBMM/reference/mbmm_run_simulation_grid.md)
  : Run a manuscript-style simulation grid
- [`mbmm_productive_infection()`](https://linhdo154.github.io/MBMM/reference/mbmm_productive_infection.md)
  : Fit a productive-infection MBMM analysis

## Model fit objects

- [`summary(`*`<mbmm_fit>`*`)`](https://linhdo154.github.io/MBMM/reference/summary.mbmm_fit.md)
  : Summarize an mbmm_fit object
- [`print(`*`<mbmm_fit>`*`)`](https://linhdo154.github.io/MBMM/reference/print.mbmm_fit.md)
  : Print an mbmm_fit object
- [`plot(`*`<mbmm_fit>`*`)`](https://linhdo154.github.io/MBMM/reference/plot.mbmm_fit.md)
  : Plot the bank log-likelihood-ratio diagnostic for an MBMM fit

## Numerical utilities

- [`mbmm_rk4_solve()`](https://linhdo154.github.io/MBMM/reference/mbmm_rk4_solve.md)
  : Fixed-step 4th-order Runge-Kutta ODE integrator (base R, no
  dependencies)
- [`mbmm_logsumexp()`](https://linhdo154.github.io/MBMM/reference/mbmm_logsumexp.md)
  : Numerically stable log-sum-exp
- [`mbmm_logmeanexp()`](https://linhdo154.github.io/MBMM/reference/mbmm_logmeanexp.md)
  : Numerically stable log-mean-exp
- [`mbmm_logit()`](https://linhdo154.github.io/MBMM/reference/mbmm_logit.md)
  : Logit
- [`mbmm_inv_logit()`](https://linhdo154.github.io/MBMM/reference/mbmm_inv_logit.md)
  : Inverse logit
- [`mbmm_safe_log10()`](https://linhdo154.github.io/MBMM/reference/mbmm_safe_log10.md)
  : Safe log10 with a positive floor
- [`mbmm_dnorm_censored()`](https://linhdo154.github.io/MBMM/reference/mbmm_dnorm_censored.md)
  : Gaussian log-density with optional left censoring at an LOD
- [`mbmm_default_distance()`](https://linhdo154.github.io/MBMM/reference/mbmm_default_distance.md)
  : Default ABC discrepancy: censoring-aware weighted RMSE on the log10
  scale
