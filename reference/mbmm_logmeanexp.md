# Numerically stable log-mean-exp

Used to Monte Carlo marginalize a class-conditional likelihood over a
reference bank of kinetic-parameter particles or bundles:

## Usage

``` r
mbmm_logmeanexp(v)
```

## Arguments

- v:

  Numeric vector of per-particle log-likelihoods.

## Value

Scalar log-mean-exp.

## Details

\$\$ \log\left\\ \frac{1}{M} \sum\_{m=1}^{M}\exp(v_m) \right\\. \$\$

Importantly, `-Inf` values remain part of the denominator \\M\\; they
represent valid particles with zero likelihood and must not be
discarded.
