# Fit a productive-infection MBMM analysis

Fits an MBMM analysis in which the productive-infection class is
represented by a mechanistic trajectory bank and the non-establishment
class is represented by an analytic observation model below the assay
LOD.

## Usage

``` r
mbmm_productive_infection(
  obs_list,
  productive_bank,
  sigma,
  formula,
  data,
  lod_fraction = 0.5,
  productive_class = "productive",
  nonestablishment_class = "nonestablishment",
  ridge = 1e-04,
  penalty = "none",
  n_perm = 0,
  seed = NULL,
  id_col = "id",
  donor_col = "id",
  bundle_col = "bundle_id",
  log_scale = FALSE,
  donor_weighting = c("equal", "particle"),
  loo_productive = FALSE,
  verbose = TRUE,
  maxit = 4000
)
```

## Arguments

- obs_list:

  Named list of subject longitudinal observations.

- productive_bank:

  Predicted productive-infection reference bank.

- sigma:

  Measurement-error SD on the Gaussian likelihood scale.

- formula:

  One-sided gate-model formula.

- data:

  Subject-level gate covariate data.frame.

- lod_fraction:

  Expected non-establishment mean as a fraction of LOD.

- productive_class:

  Label for the productive-infection class.

- nonestablishment_class:

  Label for the non-establishment class.

- ridge:

  Nonnegative ridge stabilization parameter.

- penalty:

  Either `"ridge"` or `"none"`.

- n_perm:

  Number of optional gate permutations.

- seed:

  Optional permutation-test seed.

- id_col:

  Subject-ID column in `data`.

- donor_col:

  Donor-ID column in `productive_bank`.

- bundle_col:

  Bundle-ID column in `productive_bank`.

- log_scale:

  Whether observations, predictions, and LOD values are already on the
  log10 likelihood scale.

- donor_weighting:

  Either `"equal"` or `"particle"`.

- loo_productive:

  Whether to exclude the scored subject from the productive bank when
  that subject is also a bank donor.

- verbose:

  Print scoring progress.

- maxit:

  Maximum gate optimizer iterations.

## Value

An object of class `"mbmm_fit"`.

## Details

This structure is appropriate for infection-establishment applications
in which productive infection produces structured longitudinal kinetics,
whereas non-establishment or transient/abortive detection is represented
by a low-level analytic likelihood rather than a second ODE trajectory
bank.

The expected measurement under non-establishment is defined as a
specified fraction of the observation-specific LOD. A value of `0.5`
corresponds to an expected natural-scale measurement of LOD divided by
two.

For binary models, gate coefficients describe the log-odds of the
`productive_class` relative to the `nonestablishment_class`.

## Examples

``` r
if (FALSE) { # \dontrun{
# See vignette("mbmm-intro") for a complete productive-infection example.
} # }
```
