# Analytic non-establishment class likelihood for COVID5001-style analyses

Computes the class-conditional likelihood for a non-establishment or
nonproductive class whose expected viral load is a fixed fraction of the
observation-specific LOD. In the COVID5001 manuscript analysis, the
default mean was LOD / 2 on the natural viral-load scale.

## Usage

``` r
mbmm_noise_floor_loglik(obs, sigma, lod_fraction = 0.5, log_scale = FALSE)
```

## Arguments

- obs:

  Subject observations with columns `y`, `lod`, and optionally `target`
  and `censored`.

- sigma:

  Measurement-error SD on the log10 likelihood scale. Either a scalar or
  a target-named numeric vector.

- lod_fraction:

  Expected non-establishment mean as a fraction of the LOD. The
  manuscript default is `0.5`.

- log_scale:

  Whether `obs$y` and `obs$lod` are already on the log10 scale. When
  `FALSE`, the non-establishment mean is constructed as
  `lod_fraction * lod` on the natural scale and then log10-transformed.

## Value

A scalar log-likelihood under the analytic non-establishment model.

## Details

This class is analytic and does not require a reference trajectory bank.
