# Default ABC discrepancy: censoring-aware weighted RMSE on the log10 scale

Detected observations contribute squared error on the log10 scale.
Left-censored observations are penalized only when the predicted
trajectory exceeds the observation-specific LOD.

## Usage

``` r
mbmm_default_distance(
  y_obs,
  pred,
  censored = NULL,
  lod = NULL,
  w_detect = 5,
  w_censored = 1,
  floor = 0.01
)
```

## Arguments

- y_obs:

  Observed values on the natural scale.

- pred:

  Predicted values on the natural scale.

- censored:

  Logical or 0/1 vector indicating left-censored observations.

- lod:

  Observation-specific LOD on the natural scale.

- w_detect:

  Weight for detected-observation squared errors.

- w_censored:

  Weight for censored-observation squared errors.

- floor:

  Smallest positive value passed to log10.

## Value

A nonnegative scalar discrepancy.
