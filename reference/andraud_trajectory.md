# Reduced prime-boost antibody kinetics

A reduced form of the plasma-cell and antibody model of Andraud et al.
(2012), in which antibody is produced by a short-lived and a long-lived
source and cleared at rate \\m_A\\: \$\$\dot S = -m_s S, \qquad \dot A =
W_l + S - m_A A,\$\$ with \\S(0) = W_s\\ and \\A(0) = A_0\\. Solving
gives a sum of three exponentials, and at the boost the two fluxes
receive increments, so the trajectory is continuous there by
construction.

## Usage

``` r
andraud_trajectory(
  t,
  A0,
  C1,
  C2,
  D1,
  D2,
  ms,
  mA = log(2)/(26/30),
  t_boost = 84/30,
  floor_at = 0.01
)
```

## Arguments

- t:

  Times in months since the prime.

- A0, C1, C2, D1, D2:

  Levels on the titre scale: baseline, prime-only plateau and its
  short-lived amplitude, post-boost plateau and its short-lived
  amplitude.

- ms:

  Decay rate of the short-lived source, per month.

- mA:

  Antibody clearance rate, per month. The default corresponds to an IgG
  half-life of about 26 days.

- t_boost:

  Time of the boost, months since the prime.

- floor_at:

  Value at which non-positive predictions are floored. A sum of
  exponentials with mixed signs is not constrained positive; pass `-Inf`
  to inspect the raw solution, which is what a positivity screen needs.

## Value

A numeric vector of antibody concentration at `t`.

## Details

The three levels the classes are built on are read straight off the
model. \\A_0\\ is the initial condition; setting \\\dot A = 0\\ once the
short-lived source is spent gives the prime-only plateau \\C_1 = W_l /
m_A\\; the same argument after the boost gives \\D_1 = (W_l + \Delta
W_l)/m_A\\. The transient peak four weeks after the boost reflects
\\\Delta W_s\\ and decays, so it is not the durable level.

Holding \\W_l\\ constant is what makes the plateaus well defined. It is
reasonable over a study of this length: antibody clearance is fast, yet
antibody to vaccine antigens persists for decades, which given that
clearance rate requires continuing production.

## References

Andraud M et al. (2012) Living on three time scales: the dynamics of
plasma cell and antibody populations. *PLoS Comput Biol* 8, e1002418.

## See also

[`amp_to_flux()`](https://linhdo154.github.io/MBMM/reference/amp_to_flux.md),
[`andraud_levels()`](https://linhdo154.github.io/MBMM/reference/andraud_levels.md)

## Examples

``` r
round(andraud_trajectory(c(0, 1, 2.8, 3.7, 14), A0 = 400, C1 = 1200,
                         C2 = 0, D1 = 3000, D2 = 2500, ms = 0.3))
#> [1]  400  840 1115 2774 3086
```
