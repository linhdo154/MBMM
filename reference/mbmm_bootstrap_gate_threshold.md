# Bootstrap an MBMM gate threshold

Estimates a fitted gate threshold and its uncertainty by resampling
subjects, refitting the MBMM gate, and recalculating the threshold in
each bootstrap sample.

## Usage

``` r
mbmm_bootstrap_gate_threshold(
  logL,
  formula,
  data,
  focal,
  classes = NULL,
  hold_at = list(),
  level = 0.5,
  center = 0,
  scale = 1,
  transform = c("identity", "log10", "log"),
  observed_range = NULL,
  n_boot = 1000,
  conf = 0.95,
  ridge = 1e-04,
  penalty = "ridge",
  id_col = "id",
  seed = NULL,
  maxit = 4000,
  verbose = TRUE
)
```

## Arguments

- logL:

  Data.frame containing the class-conditional log-likelihood columns
  used by
  [`mbmm_fit_gate`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md).

- formula:

  One-sided gate-model formula.

- data:

  Subject-level covariate data.frame.

- focal:

  Name of the continuous covariate whose threshold is requested.

- classes:

  Optional latent-class order passed to
  [`mbmm_fit_gate`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md).

- hold_at:

  Named list fixing all other covariates in the gate formula.

- level:

  Gate probability at which the threshold is calculated.

- center:

  Center used to standardize the focal covariate.

- scale:

  Scale used to standardize the focal covariate.

- transform:

  Transformation applied before standardization. One of `"identity"`,
  `"log10"`, or `"log"`.

- observed_range:

  Optional two-element numeric vector giving the observed range of the
  focal covariate on the fitted model scale. By default, the range is
  calculated from `data[[focal]]`.

- n_boot:

  Number of subject-level bootstrap samples.

- conf:

  Confidence level for the percentile interval.

- ridge:

  Ridge parameter passed to
  [`mbmm_fit_gate`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md).

- penalty:

  Either `"ridge"` or `"none"`.

- id_col:

  Subject-ID column in `data`.

- seed:

  Optional bootstrap seed.

- maxit:

  Maximum gate optimizer iterations.

- verbose:

  Print bootstrap progress.

## Value

An object of class `"mbmm_threshold_summary"` containing the point
estimate, percentile confidence interval, bootstrap estimates, bootstrap
success rate, and extrapolation percentage.

## Details

The function reports thresholds both on the covariate scale used in the
gate model and, optionally, on the original biomarker scale.

A bootstrap estimate is classified as extrapolated when its model-scale
threshold lies outside the observed range of the focal gate covariate.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- mbmm_bootstrap_gate_threshold(
  logL = logL,
  formula = ~ biomarker_z,
  data = biomarker_data,
  focal = "biomarker_z",
  classes = c("productive", "nonestablishment"),
  center = mean(log10(biomarker_data$biomarker)),
  scale = sd(log10(biomarker_data$biomarker)),
  transform = "log10",
  n_boot = 1000,
  seed = 1
)

result
} # }
```
