# Fit a bank-based MBMM analysis

Convenience wrapper that first scores every subject against
class-specific mechanistic reference banks and then fits the
covariate-linked MBMM gate.

## Usage

``` r
mbmm(
  obs_list,
  bank_pred_list,
  sigma,
  formula,
  data,
  classes = NULL,
  ridge = 1e-04,
  penalty = "ridge",
  n_perm = 0,
  seed = NULL,
  id_col = "id",
  donor_col = "id",
  bundle_col = "bundle_id",
  log_scale = TRUE,
  donor_weighting = c("equal", "particle"),
  loo = TRUE,
  verbose = TRUE,
  maxit = 4000
)
```

## Arguments

- obs_list:

  Named list of subject observation data.frames. Names must contain
  unique subject IDs.

- bank_pred_list:

  Named list of class-specific bank-prediction tables. Names define the
  latent-class labels.

- sigma:

  Measurement-error SD passed to
  [`mbmm_class_loglik`](https://linhdo154.github.io/MBMM/reference/mbmm_class_loglik.md).

- formula:

  One-sided gate-model formula, such as `~ baseline_antibody` or
  `~ arm`.

- data:

  Covariate data.frame containing one row per subject.

- classes:

  Optional character vector fixing the latent-class order. Defaults to
  `names(bank_pred_list)`.

- ridge:

  Nonnegative ridge penalty used when `penalty = "ridge"`.

- penalty:

  Either `"ridge"` or `"none"`.

- n_perm:

  Number of optional gate permutations.

- seed:

  Optional permutation-test seed.

- id_col:

  Name of the subject-ID column in `data`.

- donor_col:

  Name of the donor-ID column in the bank tables.

- bundle_col:

  Name of the bundle-ID column in the bank tables.

- log_scale:

  Whether observations, predictions, and LOD values are already on the
  log10 likelihood scale.

- donor_weighting:

  Either `"equal"` or `"particle"`. Equal donor weighting is the
  manuscript-oriented default.

- loo:

  Logical scalar or class-named logical vector controlling
  leave-one-subject/donor-out scoring.

- verbose:

  Print scoring progress.

- maxit:

  Maximum gate optimizer iterations.

## Value

An object of class `"mbmm_fit"` containing the scored class-conditional
likelihoods and fitted gate.

## Details

This wrapper is appropriate when every latent class is represented by a
mechanistic trajectory bank, such as a boost-responsive versus
weak/stable phenotype analysis.

For productive-infection analyses in which the non-establishment class
uses an analytic likelihood rather than a second trajectory bank, use
[`mbmm_productive_infection`](https://linhdo154.github.io/MBMM/reference/mbmm_productive_infection.md).

For binary models, gate coefficients describe the log-odds of the first
class in `classes` relative to the second class.

## Examples

``` r
obs_list <- list(
  S1 = data.frame(time = 0:1, y = c(1.0, 1.1)),
  S2 = data.frame(time = 0:1, y = c(1.2, 0.9)),
  S3 = data.frame(time = 0:1, y = c(1.1, 1.0)),
  S4 = data.frame(time = 0:1, y = c(4.8, 5.1)),
  S5 = data.frame(time = 0:1, y = c(5.2, 4.9)),
  S6 = data.frame(time = 0:1, y = c(5.0, 5.3))
)

bank_low <- data.frame(
  id = "donorA",
  bundle_id = 1,
  time = 0:1,
  pred = c(1, 1)
)

bank_high <- data.frame(
  id = "donorB",
  bundle_id = 1,
  time = 0:1,
  pred = c(5, 5)
)

covariates <- data.frame(
  id = paste0("S", 1:6),
  x = c(-1.2, -0.8, -0.4, 0.4, 0.8, 1.2)
)

fit <- mbmm(
  obs_list = obs_list,
  bank_pred_list = list(
    low = bank_low,
    high = bank_high
  ),
  sigma = 0.3,
  formula = ~x,
  data = covariates,
  classes = c("low", "high"),
  id_col = "id",
  loo = FALSE,
  verbose = FALSE
)

fit
#> MBMM fit
#>   classes: low vs high 
#>   scoring: bank-versus-bank 
#>   n subjects scored: 6 
#> MBMM gate model
#>   classes: low vs high 
#>   method: direct-marginal-optimization 
#>   penalty: ridge 
#>   ridge: 1e-04 
#>   n subjects: 6 
#>   optimizer convergence: full = 0 , null = 0 
#>   unpenalized logLik: full = 1.97207 , null = -2.18292 
#>   LRT = 8.31 on df = 1  -> p = 0.003943 
#> 
#> Gate coefficients:
#>         term      estimate   odds_ratio
#>  (Intercept) -1.289454e-11 1.000000e+00
#>            x -1.560373e+01 1.672573e-07
```
