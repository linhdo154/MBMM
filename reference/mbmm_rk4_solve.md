# Fixed-step 4th-order Runge-Kutta ODE integrator (base R, no dependencies)

A small, dependency-free ODE solver intended for building the
`simulate_fn` ingredient MBMM needs (see
[`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md)).
For stiff systems or very fine control, use deSolve instead – this
function exists so the package (and its examples) do not require any
package beyond base R.

## Usage

``` r
mbmm_rk4_solve(f, y0, times, parms, n_substeps = 10)
```

## Arguments

- f:

  function(t, y, parms) returning the derivative dy/dt as a numeric
  vector the same length as y

- y0:

  numeric vector of initial states

- times:

  numeric vector of output times (need not be evenly spaced); the first
  element is treated as t0 with state y0

- parms:

  a list or named vector passed through to `f`

- n_substeps:

  number of RK4 substeps to take between each pair of consecutive
  `times` (higher = more accurate, slower)

## Value

a matrix with one row per element of `times` and one column per state
variable (plus a leading `time` column)
