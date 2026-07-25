# Numerically stable log-sum-exp

Computes `log(sum(exp(v)))` without numerical overflow or underflow.
Negative infinity values are retained as valid zero-likelihood
contributions. Missing values are not silently removed.

## Usage

``` r
mbmm_logsumexp(v)
```

## Arguments

- v:

  Numeric vector of log-scale values.

## Value

A scalar equal to `log(sum(exp(v)))`.
