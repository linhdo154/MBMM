# Leave-one-out banks

Refits the bank once per scored participant, each time excluding that
participant's own measurements. This restores the honesty condition when
no external donor sample exists, at the cost of one fit per participant.

## Usage

``` r
bank_loo(loglik, lambda = 1, nu = NULL, ...)
```

## Arguments

- loglik:

  A participant-by-candidate log-likelihood matrix.

- lambda:

  Entropy penalty passed to
  [`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md).

- nu:

  Optional per-participant candidate retention count. When supplied,
  [`screen_candidates()`](https://linhdo154.github.io/MBMM/reference/screen_candidates.md)
  is rerun inside every leave-one-out fold using only that fold's
  donors. This is required when support screening is data-dependent;
  screening once on the full cohort before calling `bank_loo()` would
  leak the held-out participant into its bank.

- ...:

  Further arguments to
  [`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md).

## Value

A list with `weights`, a participant-by-candidate matrix whose row `i`
holds the bank fitted without participant `i`, and the corresponding
`logLik`, a vector of each participant's log class likelihood under
their own leave-one-out bank, and `support_size`, the number of
candidates retained in each fold.

## Details

The columns of `loglik` must come from a candidate family fixed
independently of the participant being scored. If the family is reduced
using participant measurements, pass `nu` here so screening is repeated
within each fold. Positivity or other screens that depend only on the
mechanistic design, and not on participant measurements, may be applied
before this function.

## See also

[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md)
