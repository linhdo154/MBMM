# Participant-by-candidate log-likelihood matrix

Scores every participant against every candidate trajectory, giving the
matrix \\\log p_k(Y_i \mid \theta_m)\\ that
[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md)
and
[`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md)
consume. This is the only place the mechanistic model enters:
`trajectory` can be any function of time and a parameter row.

## Usage

``` r
trajectory_loglik(
  obs,
  candidates,
  trajectory,
  sigma,
  lags = NULL,
  lag_prior = NULL
)
```

## Arguments

- obs:

  A data frame with one row per measurement and columns `id`, `time`,
  `value`, and optionally `censored` and `lod`.

- candidates:

  A data frame or matrix of candidate parameter values, one row per
  candidate.

- trajectory:

  A function `(time, theta)` returning the trajectory value at `time`
  for the single candidate `theta` (one row of `candidates`).

- sigma:

  Residual standard deviation on the \\\log\_{10}\\ scale.

- lags:

  Optional numeric vector of candidate time shifts, in the units of
  `obs$time`. When the clock on which measurements are recorded is not
  the biological clock – days since a reported symptom onset rather than
  days since infection – the offset is unknown and must be integrated
  out rather than fixed. Each observation time is evaluated at
  `time + lag`.

- lag_prior:

  Prior weights over `lags`, normalised internally. Defaults to uniform.

## Value

A numeric matrix with one row per participant (named by `id`) and one
column per candidate, with any lag already marginalised out.

## See also

[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md),
[`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md),
[`background_loglik()`](https://linhdo154.github.io/MBMM/reference/background_loglik.md)

## Examples

``` r
obs <- data.frame(id = rep(1:2, each = 3), time = rep(c(1, 3, 6), 2),
                  value = c(100, 60, 20, 500, 300, 90))
cand <- data.frame(rate = c(0.1, 0.3), peak = c(200, 600))
decay <- function(t, th) th$peak * exp(-th$rate * t)
trajectory_loglik(obs, cand, decay, sigma = 0.3)
#>         [,1]       [,2]
#> 1 -3.4072977 -6.2245675
#> 2 -0.7899403  0.7858733

## with an unknown infection time, integrated over a one-day grid
trajectory_loglik(obs, cand, decay, sigma = 0.3, lags = -2:2)
#>         [,1]      [,2]
#> 1 -3.0718806 -3.173365
#> 2 -0.7558832  0.337751
```
