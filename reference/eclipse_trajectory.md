# Eclipse-phase target-cell-limited viral kinetics

The within-host model used in the SARS-CoV-2 application: \$\$\dot U =
-\beta U V, \quad \dot I_1 = \beta U V - \kappa I_1, \quad \dot I_2 =
\kappa I_1 - \delta I_2, \quad \dot V = p I_2 - cV,\$\$ with \\U(0) =
T_0\\, \\I_1(0) = I_2(0) = 0\\ and \\V(0) = V_0\\. The eclipse
compartment \\I_1\\ is what separates this from the basic
target-cell-limited model: without it a newly infected cell produces
virus immediately, which the observed rise is too slow to support.

## Usage

``` r
eclipse_trajectory(
  times,
  beta,
  delta,
  p,
  V0,
  kappa = 4,
  c_clear = 10,
  T0 = 1e+07,
  grid = seq(0, max(times) + 1, by = 0.1)
)
```

## Arguments

- times:

  Numeric vector of times (days) at which to return viral load.

- beta, delta, p, V0:

  Free parameters: infection rate constant, infected-cell loss rate,
  virion production rate, inoculum.

- kappa, c_clear, T0:

  Fixed constants. Defaults are the values adopted by Ke et al. (2021)
  for the upper respiratory tract. `T0` is a scale convention rather
  than a literature value: it enters only through the products \\\beta
  T_0\\ and \\p T_0\\, so rescaling it and compensating in `beta` and
  `p` leaves the trajectory unchanged.

- grid:

  Internal solver grid.

## Value

A numeric vector of viral load at `times`.

## References

Baccam P et al. (2006) Kinetics of influenza A virus infection in
humans. *J Virol* 80, 7590-7599.

Ke R et al. (2021) In vivo kinetics of SARS-CoV-2 infection and its
relationship with a person's infectiousness. *PNAS* 118, e2111477118.

## See also

[`eclipse_design()`](https://linhdo154.github.io/MBMM/reference/eclipse_design.md)

## Examples

``` r
if (requireNamespace("deSolve", quietly = TRUE)) {
  v <- eclipse_trajectory(0:20, beta = 3e-8, delta = 1.7, p = 45, V0 = 10)
  round(log10(v[c(1, 6, 11)]), 2)
}
#> [1]  1.00 -0.52 -1.02
```
