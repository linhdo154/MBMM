# Reject candidates whose trajectory is not positive

A sum of exponentials with mixed signs is not constrained to be
positive, so a candidate can dip below zero between observations. The
likelihood sees only the observation times, so such an excursion costs
nothing there and the candidate survives into the fitted population,
where it appears in posterior summaries as a dive to the floor and back.

## Usage

``` r
screen_positive(
  candidates,
  trajectory,
  grid = seq(0, 60, by = 1),
  floor = 0.01
)
```

## Arguments

- candidates:

  A data frame of candidate parameter values.

- trajectory:

  A function `(time, theta)` as in
  [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md).
  It must be called in a way that returns the raw solution, not one
  already floored at a positive value.

- grid:

  Dense time grid on which positivity is checked.

- floor:

  Smallest value treated as positive.

## Value

A list with the logical vector `keep`, the screened `candidates`, and
the `rejected` count.

## Details

The grid must be dense enough to catch an excursion *between* visits;
evaluating on the observation times defeats the purpose.

## See also

[`screen_candidates()`](https://linhdo154.github.io/MBMM/reference/screen_candidates.md)

## Examples

``` r
cand <- data.frame(a = c(1, 1), b = c(0.5, 5))
f <- function(t, th) th$a - th$b * t          # the second goes negative
screen_positive(cand, f, grid = seq(0, 1, by = 0.05))$keep
#> [1]  TRUE FALSE
```
