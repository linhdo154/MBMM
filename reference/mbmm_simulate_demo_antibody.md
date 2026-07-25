# Simulate a fully synthetic baseline antibody measurement

The public demonstration uses a bounded two-component Gaussian mixture
on the natural-log scale. This produces a broad, right-skewed antibody
distribution on the natural scale without using empirical COVID5001
measurements.

## Usage

``` r
mbmm_simulate_demo_antibody(
  n,
  high_fraction = 0.25,
  low_mean = 1.25,
  high_mean = 2.35,
  low_sd = 0.55,
  high_sd = 0.4,
  clip = c(0, 3.3),
  calibration_mean = NULL,
  calibration_sd = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of subjects.

- high_fraction:

  Probability of belonging to the higher-response component.

- low_mean, high_mean:

  Component means on the natural-log scale.

- low_sd, high_sd:

  Component standard deviations on the natural-log scale.

- clip:

  Lower and upper bounds applied on the natural-log scale.

- calibration_mean:

  Optional fixed mean used to standardize the biomarker.

- calibration_sd:

  Optional fixed standard deviation used to standardize it.

- seed:

  Optional random seed.

## Value

A data frame containing natural-scale, log-scale, and standardized
antibody values.
