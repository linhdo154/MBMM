# Create the default reduced manuscript-style scenario grid

This grid is deliberately tiny and is intended for a vignette
demonstration, not for estimating statistical power or type-I error.

## Usage

``` r
mbmm_demo_scenario_grid(
  replicates = 2L,
  n_subjects = 50L,
  prevalence = 0.3,
  beta_grid = c(0, -1.5),
  p_abortive_grid = c(0, 0.4),
  base_seed = 20260724L
)
```

## Arguments

- replicates:

  Number of demonstration replicates per setting.

- n_subjects:

  Number of subjects per generated dataset.

- prevalence:

  Target marginal productive-infection prevalence.

- beta_grid:

  True biomarker effects.

- p_abortive_grid:

  Abortive-pattern probabilities.

- base_seed:

  First seed in the grid.

## Value

A data frame.
