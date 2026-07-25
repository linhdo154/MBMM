# Predict mechanistic bank trajectories at required observation times

Applies a user-supplied mechanistic simulator to every particle or
multitarget bundle row in a reference bank.

## Usage

``` r
mbmm_predict_bank(
  bank,
  simulate_fn,
  times = NULL,
  times_by_target = NULL,
  donor_col = "id",
  bundle_col = "bundle_id",
  failure = c("error", "drop")
)
```

## Arguments

- bank:

  Reference-bank data.frame.

- simulate_fn:

  Function with signature `simulate_fn(theta_row, times)` returning one
  numeric prediction per requested time.

- times:

  Numeric vector of evaluation times for single-target banks.

- times_by_target:

  Named list of time vectors for multitarget banks.

- donor_col:

  Name of the donor ID column.

- bundle_col:

  Name of the bundle ID column.

- failure:

  Either `"error"` or `"drop"`. With `"error"`, any failed or invalid
  simulation stops. With `"drop"`, failed bank rows are omitted and
  reported in attributes.

## Value

Long-format bank predictions containing donor ID, bundle ID, optional
target, time, and prediction.
