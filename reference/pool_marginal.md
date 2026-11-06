# Integrate a nuisance pool out of a candidate likelihood

A mechanistic design often splits into coordinates the question is about
and coordinates it is not. In the prime-boost application the class
coordinate \\g\\ and the decay rate \\m_s\\ are shared and carry the
question, while the per-assay levels \\(A_0, C_1)\\ and the spike ratios
are nuisance: they differ between antigen targets and nothing is being
asked about them.

## Usage

``` r
pool_marginal(loglik, shared, weights = NULL)
```

## Arguments

- loglik:

  A participant-by-candidate log-likelihood matrix, as returned by
  [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md).

- shared:

  A vector of length `ncol(loglik)` giving the shared point each
  candidate belongs to. A factor fixes the output column order; anything
  else is converted with
  [`factor()`](https://rdrr.io/r/base/factor.html).

- weights:

  Optional pool weights, length `ncol(loglik)`, normalised within each
  shared point. The default is an equal-weight average, which is the
  Monte Carlo estimate of the integral when the pool is drawn from the
  nuisance prior.

## Value

A matrix with one row per participant and one column per shared point,
holding \\\log \overline{p}(Y_i \mid s)\\. Column names are the shared
point labels.

## Details

Those nuisance coordinates must be **integrated out**, not competed
against. If they are left as separate bank atoms, a participant can be
matched by any candidate that happens to fit, whatever its \\g\\, and
the class coordinate stops being identified — a flat trajectory then
scores equally well at \\g\>0\\ and \\g\<0\\. Averaging instead asks the
right question: *given this shared point, how well does a typical draw
of the nuisance coordinates explain this participant?* A flat
participant scores badly at \\g\>0\\ because most of the pool predicts a
rise there.

The result is the inner integral of the two-level construction, and is
what
[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md)
should be given: the bank's support is the shared grid, not the
candidate list.

## See also

[`jitter_marginal()`](https://linhdo154.github.io/MBMM/reference/jitter_marginal.md)
for integrating a calibration offset,
[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md)
for the outer problem.

## Examples

``` r
ll <- matrix(c(-1, -2, -3, -4, -5, -6), nrow = 2)
pool_marginal(ll, shared = c("a", "a", "b"))
#>              a  b
#> [1,] -1.566219 -5
#> [2,] -2.566219 -6
```
