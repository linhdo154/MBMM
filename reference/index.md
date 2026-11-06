# Package index

## Likelihood

Turning measurements into evidence for a state

- [`obs_loglik()`](https://linhdo154.github.io/MBMM/reference/obs_loglik.md)
  : Log-likelihood of one measurement under a trajectory
- [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md)
  : Participant-by-candidate log-likelihood matrix
- [`background_loglik()`](https://linhdo154.github.io/MBMM/reference/background_loglik.md)
  : Background (class-0) log-likelihood at a fraction of the assay limit
- [`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md)
  : Class likelihood from a bank

## Building the likelihood over shared points

Integrating out what the question is not about

- [`pool_marginal()`](https://linhdo154.github.io/MBMM/reference/pool_marginal.md)
  : Integrate a nuisance pool out of a candidate likelihood
- [`jitter_marginal()`](https://linhdo154.github.io/MBMM/reference/jitter_marginal.md)
  : Integrate a calibration offset along a grid coordinate

## The kinetic population

Estimating the bank, and keeping it honest

- [`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md)
  : Estimate bank weights by entropy-penalised maximum likelihood
- [`bank_loo()`](https://linhdo154.github.io/MBMM/reference/bank_loo.md)
  : Leave-one-out banks
- [`screen_candidates()`](https://linhdo154.github.io/MBMM/reference/screen_candidates.md)
  : Screen candidates to the ones participants actually support
- [`screen_positive()`](https://linhdo154.github.io/MBMM/reference/screen_positive.md)
  : Reject candidates whose trajectory is not positive
- [`profile_sigma()`](https://linhdo154.github.io/MBMM/reference/profile_sigma.md)
  : Profile the residual scale
- [`restrict_bank()`](https://linhdo154.github.io/MBMM/reference/restrict_bank.md)
  : Restrict a bank to a region of parameter space (Route II)

## The gate

Covariate association with the latent state

- [`fit_gate()`](https://linhdo154.github.io/MBMM/reference/fit_gate.md)
  : Fit the covariate gate by EM
- [`gate_test()`](https://linhdo154.github.io/MBMM/reference/gate_test.md)
  : Test a covariate effect in the gate

## Mechanistic models

- [`eclipse_trajectory()`](https://linhdo154.github.io/MBMM/reference/eclipse_trajectory.md)
  : Eclipse-phase target-cell-limited viral kinetics
- [`eclipse_design()`](https://linhdo154.github.io/MBMM/reference/eclipse_design.md)
  : Candidate design on identifiable viral coordinates
- [`andraud_trajectory()`](https://linhdo154.github.io/MBMM/reference/andraud_trajectory.md)
  : Reduced prime-boost antibody kinetics
- [`andraud_levels()`](https://linhdo154.github.io/MBMM/reference/andraud_levels.md)
  : Class coordinate for a prime-boost question
- [`amp_to_flux()`](https://linhdo154.github.io/MBMM/reference/amp_to_flux.md)
  [`flux_to_amp()`](https://linhdo154.github.io/MBMM/reference/amp_to_flux.md)
  : Convert between amplitudes and production fluxes

## Simulation

- [`simulate_mbmm()`](https://linhdo154.github.io/MBMM/reference/simulate_mbmm.md)
  : Simulate a cohort with a latent kinetic state

## Visualisation

Seeing the kinetic population and checking simulation recovery

- [`plot_bank()`](https://linhdo154.github.io/MBMM/reference/plot_bank.md)
  : Plot a weighted mechanistic trajectory bank
- [`plot_responsibility()`](https://linhdo154.github.io/MBMM/reference/plot_responsibility.md)
  : Plot posterior responsibilities against simulated truth
- [`plot_gate()`](https://linhdo154.github.io/MBMM/reference/plot_gate.md)
  : Plot a fitted covariate gate
