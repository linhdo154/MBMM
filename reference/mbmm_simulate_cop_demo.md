# Simulate a synthetic correlate-of-protection cohort

Generates a small synthetic dataset in the spirit of the simulation
study in Sec. 3.2 of Do, McCarthy & Rogers: each subject has a baseline
(standardized) antibody value; a latent productive-infection state `Z`
is drawn from an antibody-dependent logistic gate; subjects with `Z = 1`
get a viral-load trajectory from the eclipse-phase model
([`mbmm_example_eclipse_simulate`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_simulate.md))
under subject-specific true kinetic parameters, observed with log-normal
noise and a limit of detection; subjects with `Z = 0` get a flat
near-LOD trajectory with an occasional single-point transient/abortive
"blip" (mimicking a non-productive PCR-detectable event). A subset of
the `Z = 1` subjects is given a denser visit schedule and flagged as
bank-eligible, mirroring the donor-richness threshold used for real
reference banks (Sec. 4.3).

## Usage

``` r
mbmm_simulate_cop_demo(
  n_subjects = 90,
  prevalence = 0.35,
  beta_true = -1.2,
  n_bank_donors = 18,
  sigma_obs = 0.35,
  lod = 3,
  seed = 1
)
```

## Arguments

- n_subjects:

  total number of subjects to simulate

- prevalence:

  marginal probability that a "typical" (zAb = 0) subject is truly
  infected (Z = 1); combined with `beta_true` to set the gate intercept

- beta_true:

  true log-odds effect of standardized baseline antibody on infection
  probability (negative = protective)

- n_bank_donors:

  number of Z = 1 subjects given a dense visit schedule and flagged
  `bank_eligible = TRUE`

- sigma_obs:

  log10-scale observation noise SD

- lod:

  limit of detection (natural scale, e.g. copies/mL)

- seed:

  RNG seed

## Value

a list with elements `obs_list` (named list of per-subject data.frames
with columns `time`, `y`, `censored`, `lod`), `covariates` (data.frame
with `id` and `baseline_antibody`), and `truth` (data.frame with the
true generating values, for checking recovery – never used by the
fitting functions themselves)

## Details

This function only \*generates data\*; it does not call any MBMM fitting
functions. See
[`vignette("mbmm-intro")`](https://linhdo154.github.io/MBMM/articles/mbmm-intro.md)
or
`system.file("examples/01_simulated_correlate_of_protection.R", package = "mbmm")`
for the full analysis pipeline built on top of it.
