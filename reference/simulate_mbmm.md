# Simulate a cohort with a latent kinetic state

Generates the setting the method is built for: a baseline marker shifts
the probability of a latent state, participants in that state follow a
mechanistic trajectory, and some participants *not* in it still produce
a brief detectable event. A rule that counts any detection as the state
therefore mislabels part of the cohort, which is what attenuates a
covariate effect and costs power.

## Usage

``` r
simulate_mbmm(
  n,
  beta_ab = -1,
  prevalence = 0.2,
  p_abortive = 0.4,
  times = c(0, 2, 4, 7, 10, 14, 21),
  lod = 100,
  sigma_obs = 0.5,
  n_pool = 200L,
  seed = NULL
)
```

## Arguments

- n:

  Number of participants.

- beta_ab:

  Covariate effect on the log-odds of the state. Negative means a higher
  marker lowers the probability.

- prevalence:

  Target marginal prevalence of the state. The gate intercept is solved
  to hit it against a large reference draw, so the realised prevalence
  does not drift with `n`.

- p_abortive:

  Fraction of participants *not* in the state who still produce a
  transient detectable event. This is the contamination the method is
  meant to absorb; at zero, a detection rule is already faithful and
  there is little for a latent-state model to recover.

- times:

  Visit schedule.

- lod:

  Assay limit of detection.

- sigma_obs:

  Observation noise, \\\log\_{10}\\ scale.

- n_pool:

  Size of the pre-solved trajectory population drawn from
  [`eclipse_design()`](https://linhdo154.github.io/MBMM/reference/eclipse_design.md).

- seed:

  Optional seed.

## Value

A list with `obs` (long-format measurements, with `censored` and `lod`),
`truth` (one row per participant: the latent `z`, the covariate `x`, and
the generating parameters), and the `alpha` that was solved for.

## Details

Because the state is generated rather than inferred, the simulation is
the one setting in which recovery can be checked directly. It is also
the only worked example the package can ship: the analyses in the
accompanying manuscript use access-controlled cohorts.

## See also

[`eclipse_design()`](https://linhdo154.github.io/MBMM/reference/eclipse_design.md),
[`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md)

## Examples

``` r
if (requireNamespace("deSolve", quietly = TRUE)) {
  sim <- simulate_mbmm(n = 12, beta_ab = -1, n_pool = 8, seed = 1)
  table(sim$truth$z)
}
#> 
#>  0 
#> 12 
```
