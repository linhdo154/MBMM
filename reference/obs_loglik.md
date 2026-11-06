# Log-likelihood of one measurement under a trajectory

Observations are modelled on the \\\log\_{10}\\ scale. A detected
measurement contributes a normal density about the trajectory; a
measurement below the assay limit contributes the corresponding tail
probability, so left-censoring is handled exactly rather than by
imputation.

## Usage

``` r
obs_loglik(y, mu, sigma, lod = NULL, censored = FALSE, eps = 1e-12)
```

## Arguments

- y:

  Observed value on the natural (not log) scale.

- mu:

  Trajectory value at the observation time, natural scale.

- sigma:

  Residual standard deviation on the \\\log\_{10}\\ scale.

- lod:

  Assay limit of detection, natural scale. Used only when `censored` is
  `TRUE`.

- censored:

  Logical; `TRUE` if the measurement is below `lod`.

- eps:

  Small constant guarding [`log10()`](https://rdrr.io/r/base/Log.html)
  against zero.

## Value

A numeric vector of log-likelihood contributions.

## Examples

``` r
obs_loglik(y = 1000, mu = 1200, sigma = 0.5)
#> [1] -0.2383307
obs_loglik(y = NA, mu = 50, sigma = 0.5, lod = 100, censored = TRUE)
#> [1] -0.3196092
```
