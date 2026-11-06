# Class likelihood from a bank

Integrates the trajectory likelihood over an estimated kinetic
population, \\\widehat L\_{ik} = \sum_m \widehat p\_{km} \\ p_k(Y_i \mid
\theta\_{km})\\. Returned on the log scale.

## Usage

``` r
class_loglik(loglik, weights)
```

## Arguments

- loglik:

  A participant-by-candidate log-likelihood matrix, as returned by
  [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md).

- weights:

  Bank weights summing to one, one per candidate.

## Value

A named numeric vector of log class likelihoods, one per participant.
