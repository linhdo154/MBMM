# Fit the covariate gate by EM

Combines the prior state probability \\\pi_i(\gamma) = \mathrm{expit}
(x_i^\top \gamma)\\ with the two class likelihoods to form the posterior
responsibility \$\$w_i(\gamma) = \frac{\pi_i L\_{i1}}{\pi_i L\_{i1} +
(1-\pi_i) L\_{i0}},\$\$ then refits \\\gamma\\ by logistic regression
using \\w_i\\ as a fractional outcome, iterating to convergence. A
participant whose trajectory does not discriminate the states keeps a
responsibility close to the gate probability and contributes at that
weight rather than being forced into a class.

## Usage

``` r
fit_gate(l1, l0, X, firth = FALSE, maxit = 2000L, tol = 1e-11)
```

## Arguments

- l1, l0:

  Log class likelihoods for state 1 and state 0, one per participant.
  See
  [`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md).

- X:

  Gate design matrix, one row per participant, including an intercept
  column.

- firth:

  Use the Jeffreys-penalised score of Firth (1993). Recommended when any
  fitted probability approaches zero or one, which is routine with few
  participants per covariate level.

- maxit, tol:

  Iteration controls.

## Value

An object of class `mbmm_gate`: `gamma`, the observed-data `loglik`, the
responsibilities `w`, fitted probabilities `pi`, and, when
`firth = TRUE`, the `penalised` objective \\\ell + \tfrac12 \log
\|X^\top W X\|\\.

## References

Firth D (1993). Bias reduction of maximum likelihood estimates.
*Biometrika* 80, 27-38.

## See also

[`gate_test()`](https://linhdo154.github.io/MBMM/reference/gate_test.md)

## Examples

``` r
set.seed(1)
n <- 60; x <- rnorm(n); X <- cbind(1, x)
l1 <- rnorm(n, 1); l0 <- rnorm(n)
fit_gate(l1, l0, X)$gamma
#>                       x 
#> 24.76340468 -0.03506096 
```
