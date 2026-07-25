# Example mechanistic trajectory model: eclipse-phase viral kinetics

A target-cell-limited within-host viral kinetic model with an eclipse
phase (Eq. 11-14 of Do, McCarthy & Rogers / Baccam et al. 2006),
provided as a ready-to-use example of the `simulate_fn` ingredient MBMM
needs (see
[`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md)).
Integration uses deSolve (`lsoda`), because this system is numerically
stiff enough that a naive fixed-step integrator (like
[`mbmm_rk4_solve`](https://linhdo154.github.io/MBMM/reference/mbmm_rk4_solve.md))
can become unstable at the eclipse/clearance rates used here; deSolve is
a package dependency for exactly this reason. Fixed constants (`k`, `c`)
and prior locations below are illustrative, not pathogen-specific –
re-derive them from the literature for a real analysis (as Sec. 4.2/App.
A.1 of the manuscript do from UKC19 human challenge data).

## Usage

``` r
mbmm_example_eclipse_simulate(theta_row, times, T0 = 1e+07, k = 4, c = 10)
```

## Arguments

- theta_row:

  a one-row data.frame/list with elements `beta`, `delta`, `p`, `V0`,
  `offset`

- times:

  observed times (e.g. "day post enrollment"), any non-negative numeric
  vector

- T0:

  fixed initial target-cell count (default 1e7)

- k:

  fixed eclipse-to-productive transition rate, 1/day (default 4)

- c:

  fixed free-virus clearance rate, 1/day (default 10)

## Value

numeric vector of predicted viral load (natural scale, e.g. copies/mL)
at each element of `times`

## Details

Free parameters (columns expected in `theta_row`): `beta` (infection
rate), `delta` (clearance rate of infected cells), `p` (viral production
rate), `V0` (initial viral burden), and `offset` (days between the
model's internal time origin and the cohort's observed day-0, allowing
unobserved infection timing to vary across subjects – see Sec. 4.1 of
the manuscript).
