# Gaussian log-density with optional left censoring at an LOD

Detected observations contribute a Gaussian log-density on the
measurement scale. Left-censored observations contribute the Gaussian
probability of falling at or below the LOD.

## Usage

``` r
mbmm_dnorm_censored(y, mean, sd, censored = FALSE, lod = NULL)
```

## Arguments

- y:

  Observed value on the same scale as `mean`. Values of `y` are ignored
  where `censored` is TRUE.

- mean:

  Predicted mean on the likelihood scale.

- sd:

  Strictly positive measurement-error standard deviation.

- censored:

  Logical or 0/1 vector indicating left-censored observations.

- lod:

  LOD on the same scale as `y` and `mean`. Required for every censored
  observation.

## Value

Numeric vector of pointwise log-likelihood contributions.

## Details

All arguments must already be on the same likelihood scale, normally the
log10 scale used by MBMM.
