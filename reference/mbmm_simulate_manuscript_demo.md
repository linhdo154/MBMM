# Simulate a reduced manuscript-style MBMM dataset

This public simulation preserves the main manuscript structure while
replacing private empirical biomarker and visit-pattern inputs with
fully synthetic analogues.

## Usage

``` r
mbmm_simulate_manuscript_demo(
  n_subjects = 50,
  prevalence = 0.3,
  beta_true = -1.5,
  p_abortive = 0.4,
  n_bank_donors = 8L,
  productive_parameter_bank = NULL,
  visit_schedules = mbmm_demo_visit_schedules(),
  sigma_obs = 0.5,
  lod = 100,
  productive_min_detect = 1L,
  productive_max_tries = 50L,
  seed = NULL
)
```

## Arguments

- n_subjects:

  Number of subjects.

- prevalence:

  Target marginal productive-infection prevalence for the realized
  synthetic biomarker distribution.

- beta_true:

  True standardized biomarker effect in the logistic gate. Negative
  values represent protection.

- p_abortive:

  Probability that a nonproductive subject receives an
  abortive/transient pattern rather than a purely near-background
  trajectory.

- n_bank_donors:

  Number of productive subjects assigned the complete visit schedule and
  flagged as bank eligible.

- productive_parameter_bank:

  Optional data frame of productive ODE parameters. When omitted, a
  synthetic public bank is generated from
  [`mbmm_example_eclipse_prior`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_prior.md).

- visit_schedules:

  Named list of public visit schedules.

- sigma_obs:

  Observation-noise SD on the log10 scale.

- lod:

  Natural-scale assay limit of detection.

- productive_min_detect:

  Minimum number of detections requested for a productive trajectory.

- productive_max_tries:

  Maximum attempts used to obtain that minimum.

- seed:

  Optional random seed.

## Value

A list with:

- `subject_data`:

  Subject-level covariates and truth.

- `truth`:

  Alias of `subject_data`.

- `observations`:

  Long-format observations.

- `obs_list`:

  Named per-subject observation data frames.

- `covariates`:

  Data frame containing `id` and `baseline_antibody`.

- `productive_parameter_bank`:

  Parameter bank used for productive trajectory generation.

- `settings`:

  Generating settings.

## Details

Productive subjects are generated from
[`mbmm_example_eclipse_simulate`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_simulate.md)
using parameter draws from
[`mbmm_example_eclipse_prior`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_prior.md)
unless a custom productive parameter bank is supplied.
