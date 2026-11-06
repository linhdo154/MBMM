# Candidate design on identifiable viral coordinates

Draws candidate trajectories on the coordinates the data determine –
doubling time, peak capacity, clearance rate and inoculum – and derives
the ODE rate constants from them. A box placed directly on \\(\beta,
\delta, p)\\ forces \\\log\_{10} K \ge 5.65\\ for every candidate, so an
established infection cannot peak low; the fitted population then
acquires a spurious non-establishing family in its place.

## Usage

``` r
eclipse_design(
  n,
  log10_K = c(3, 9.5),
  Td = c(0.3, 3),
  delta = c(0.05, 5),
  log10_V0 = c(0, 3),
  kappa = 4,
  c_clear = 10,
  T0 = 1e+07,
  seed = NULL
)
```

## Arguments

- n:

  Number of candidates.

- log10_K:

  Range of the peak capacity \\K = p T_0 / c\\, log10 scale.

- Td:

  Range of the doubling time in days, sampled log-uniformly.

- delta:

  Range of the infected-cell loss rate, per day.

- log10_V0:

  Range of the inoculum, log10 scale.

- kappa, c_clear, T0:

  Fixed constants, as in
  [`eclipse_trajectory()`](https://linhdo154.github.io/MBMM/reference/eclipse_trajectory.md).

- seed:

  Optional seed, so a design can be reproduced exactly.

## Value

A data frame of candidates with the identifiable coordinates and the
derived `beta` and `p`.

## Details

Because the doubling time is sampled positive, every candidate has \\r
\> 0\\ and hence \\R_0 \> 1\\: the class-1 support consists of
establishing kinetics **by construction** rather than by filtering
fitted output. That is a property of the design, not a finding.

## Examples

``` r
head(eclipse_design(5, seed = 1))
#>          Td            K     delta        V0         beta            p kappa
#> 1 2.3741652     53187.10 1.0695741  31.12417 2.826918e-05   0.05318710     4
#> 2 2.6411710    262304.26 0.9239559 142.18557 4.946214e-06   0.26230426     4
#> 3 1.3737859   5291110.87 3.4507631 945.62356 8.843079e-07   5.29111087     4
#> 4 1.2771306 800480275.42 1.9513134  13.80720 3.730490e-09 800.48027542     4
#> 5 0.3458657     20461.27 3.8607150 214.94912 5.164616e-04   0.02046127     4
#>   c_clear    T0
#> 1      10 1e+07
#> 2      10 1e+07
#> 3      10 1e+07
#> 4      10 1e+07
#> 5      10 1e+07
```
