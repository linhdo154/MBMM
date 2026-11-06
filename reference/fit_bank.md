# Estimate bank weights by entropy-penalised maximum likelihood

Given candidate trajectories and the log-likelihood of every participant
under each, estimates how often each candidate occurs in the population.
The objective is the marginal log-likelihood of the donor measurements
plus an entropy penalty, \$\$J\_\lambda(p) = \sum\_{j} \log \sum_m p_m
\\ p(Y_j \mid \theta_m) + \lambda H(p), \qquad H(p) = -\sum_m p_m \log
p_m .\$\$

## Usage

``` r
fit_bank(loglik, lambda = 1, start = NULL, maxit = 5000L, tol = 1e-08)
```

## Arguments

- loglik:

  A donor-by-candidate log-likelihood matrix, as returned by
  [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md).
  Donors must not include any participant who will later be scored
  against the resulting bank (see Details).

- lambda:

  Entropy penalty, non-negative.

- start:

  Optional starting weights; defaults to uniform.

- maxit, tol:

  Optimiser controls.

## Value

An object of class `mbmm_bank`: a list with `weights`, the
donor-by-candidate posterior `responsibility`, the attained `objective`,
`lambda`, and convergence information.

## Details

With `lambda = 0` this is the nonparametric maximum likelihood estimator
on a fixed support, whose solution is extremely sparse and unstable when
many near-identical trajectories compete for the same donors. The
entropy term spreads weight across such near-duplicates. For
`lambda > 0`, the entropy term is strictly concave and the maximum is
unique in the relative interior of the simplex. At `lambda = 0`, a
maximum need not be unique.

The objective is a **sum** over donors, not a mean. Raising `lambda`
with a fixed donor set therefore strengthens the penalty relative to a
fixed amount of data, but holding `lambda` fixed as donors accumulate
lets the penalty become negligible.

**Honesty.** The bank used to score a participant must be built without
that participant's measurements. Pass a disjoint donor set, or use
[`bank_loo()`](https://linhdo154.github.io/MBMM/reference/bank_loo.md)
to refit leaving each scored participant out in turn.

## See also

[`bank_loo()`](https://linhdo154.github.io/MBMM/reference/bank_loo.md),
[`restrict_bank()`](https://linhdo154.github.io/MBMM/reference/restrict_bank.md),
[`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md)

## Examples

``` r
set.seed(1)
LL <- matrix(rnorm(20 * 6), 20, 6)
b <- fit_bank(LL, lambda = 1)
round(b$weights, 3)
#> [1] 0.336 0.030 0.058 0.203 0.268 0.105
```
