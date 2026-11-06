# Test a covariate effect in the gate

Compares a full gate against a reduced gate in which the effect is
removed. With `firth = FALSE` the statistic is the usual likelihood
ratio \$\$T = 2\\\hat\ell(\hat\gamma\_{\mathrm{full}}) -
\hat\ell(\hat\gamma\_{\mathrm{red}})\\.\$\$

## Usage

``` r
gate_test(l1, l0, X, X0 = NULL, firth = FALSE, B = 0L, seed = NULL)
```

## Arguments

- l1, l0:

  Log class likelihoods, one per participant.

- X:

  Full gate design matrix, including an intercept.

- X0:

  Reduced design matrix. Defaults to the intercept alone.

- firth:

  Use the Firth-penalised gate and \\T\_{\mathrm{PO}}\\.

- B:

  Number of permutations. `0` skips the permutation null. The class
  likelihoods are held fixed while the covariate rows are shuffled,
  which is valid because they do not use the covariate.

- seed:

  Optional seed for the permutation draw.

## Value

An object of class `mbmm_test` with the observed `statistic`, its `df`,
the permutation `p_perm` and null distribution when `B > 0`, the
chi-squared `p_chisq` when `firth = FALSE`, and both fitted gates.

## Details

With `firth = TRUE` it is instead twice the difference of **penalised
objectives**, \$\$T\_{\mathrm{PO}} =
2\\\tilde\ell_X(\hat\gamma\_{\mathrm{full}}) -
\tilde\ell\_{X_0}(\hat\gamma\_{\mathrm{red}})\\, \qquad
\tilde\ell_X(\gamma) = \ell(\gamma) + \tfrac12 \log \|X^\top W(\gamma)
X\| .\$\$ Because the two penalty determinants have different
dimensions, \\T\_{\mathrm{PO}}\\ is a penalised-objective difference and
**not** a likelihood ratio. No chi-squared reference is reported for it;
use the permutation null, which is valid for any fixed statistic.

## See also

[`fit_gate()`](https://linhdo154.github.io/MBMM/reference/fit_gate.md)

## Examples

``` r
set.seed(1)
n <- 80; x <- rnorm(n); X <- cbind(1, x)
l1 <- rnorm(n, 0.6 * x); l0 <- rnorm(n)
gate_test(l1, l0, X, B = 99, seed = 1)
#> <mbmm_test>
#>   T (likelihood ratio) = 0.382 on 1 df
#>   chi-squared p = 0.5365
#>   permutation p = 0.51000  (B = 99, smallest attainable 0.01000)
```
