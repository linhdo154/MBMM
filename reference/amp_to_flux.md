# Convert between amplitudes and production fluxes

The implementation samples the levels the data speak to – \\(A_0, C_1,
C_2, D_1, D_2)\\, all in titre units – and derives the fluxes. The
fluxes are not in titre units: their scale depends on \\m_A\\ and
\\m_s\\, so a box on them induces an uncontrolled range of attainable
curves. The map is a bijection for fixed \\(m_A, m_s)\\, so this is a
reparameterisation and not a different model.

## Usage

``` r
amp_to_flux(A0, C1, C2, D1, D2, ms, mA = log(2)/(26/30), t_boost = 84/30)

flux_to_amp(A0, Ws, Wl, dWs, dWl, ms, mA = log(2)/(26/30), t_boost = 84/30)
```

## Arguments

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

- Ws, Wl, dWs, dWl:

  Production fluxes and their increments at the boost.

## Value

A list of the converted parameters.

## Examples

``` r
f <- amp_to_flux(400, 1200, 0, 3000, 2500, ms = 0.3)
str(flux_to_amp(f$A0, f$Ws, f$Wl, f$dWs, f$dWl, f$ms))
#> List of 6
#>  $ A0: num 400
#>  $ C1: num 1200
#>  $ C2: num 0
#>  $ D1: num 3000
#>  $ D2: num 2500
#>  $ ms: num 0.3
```
