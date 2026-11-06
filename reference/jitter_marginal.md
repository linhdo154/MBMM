# Integrate a calibration offset along a grid coordinate

Several assays can measure the same underlying gain on scales that do
not quite agree. Rather than give each assay its own class coordinate,
which would stop them speaking to the same question, the shared
coordinate is read with an assay-specific offset: assay \\t\\ sees \\g +
w\delta_t\\, and the offset is integrated out with equal weight.

## Usage

``` r
jitter_marginal(E, n_grid, keep, offsets, weights = NULL)
```

## Arguments

- E:

  A participant-by-column matrix whose columns are the grid coordinate
  (running fastest) crossed with any number of other coordinates.

- n_grid:

  The length of the grid coordinate, which must divide `ncol(E)`.

- keep:

  Integer indices into `seq_len(n_grid)` to retain: the interior
  sub-range on which the result is reported.

- offsets:

  Integer offsets, in grid steps, to average over. Every
  `keep + offsets` must lie inside the grid, which is checked.

- weights:

  Optional weights on `offsets`; equal by default.

## Value

A matrix with `length(keep)` times the other coordinates' columns, the
grid coordinate still running fastest.

## Details

This requires the likelihood on a *wider* grid than the one kept, so
that every shifted read stays inside it. That is why the design has an
outer grid and an interior shared sub-range: `keep` indexes the latter
within the former.

## See also

[`pool_marginal()`](https://linhdo154.github.io/MBMM/reference/pool_marginal.md).

## Examples

``` r
E <- matrix(log(c(1, 2, 3, 4, 5)), nrow = 1)
jitter_marginal(E, n_grid = 5, keep = 2:4, offsets = c(-1, 0, 1))
#>           [,1]     [,2]     [,3]
#> [1,] 0.6931472 1.098612 1.386294
```
