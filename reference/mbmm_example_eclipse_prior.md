# Example ABC prior sampler for the eclipse-phase model

Draws candidate parameters for
[`mbmm_example_eclipse_simulate`](https://linhdo154.github.io/MBMM/reference/mbmm_example_eclipse_simulate.md)
from independent priors calibrated to give a growth/peak/decline
trajectory over roughly 3 weeks (basic reproduction number \\R_0 = \beta
T_0 p / (\delta c)\\ centered in the 5-15 range at `T0 = 1e7`,
`c = 10`). Intended as a template: real analyses should replace these
with priors informed by the literature or an external challenge-study
posterior, as in Sec. 4.2/App. A.1 of the manuscript.

## Usage

``` r
mbmm_example_eclipse_prior(n)
```

## Arguments

- n:

  number of draws

## Value

a data.frame with columns `beta`, `delta`, `p`, `V0`, `offset`
