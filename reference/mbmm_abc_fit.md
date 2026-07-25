# Fit an ABC rejection-sampling posterior for one subject

Draws kinetic parameter vectors from the prior, simulates the
mechanistic trajectory, evaluates a censoring-aware discrepancy, and
retains the lowest-discrepancy candidates.

## Usage

``` r
mbmm_abc_fit(
  id,
  obs,
  simulate_fn,
  prior_sample_fn,
  n_draws = 10000,
  n_keep = NULL,
  accept_frac = 0.025,
  distance_fn = NULL,
  seed = NULL,
  failure = c("discard", "error")
)
```

## Arguments

- id:

  Subject identifier.

- obs:

  Subject observation data.frame containing `time` and `y`, and
  optionally `censored` and `lod`.

- simulate_fn:

  Function with signature `simulate_fn(theta_row, times)`.

- prior_sample_fn:

  Function with signature `prior_sample_fn(n)` returning a data.frame.

- n_draws:

  Number of prior draws.

- n_keep:

  Optional number of particles to retain.

- accept_frac:

  Optional fraction of valid draws to retain. Used when `n_keep = NULL`.

- distance_fn:

  Optional custom function with signature
  `distance_fn(theta_row, obs, pred)`.

- seed:

  Optional random seed.

- failure:

  Either `"discard"` or `"error"`. With `"discard"`, failed simulations
  receive infinite discrepancy.

## Value

A data.frame of accepted particles with `distance` and `id` columns. ABC
diagnostic information is stored as attributes.

## Details

The number retained can be specified directly using `n_keep`, or by an
acceptance fraction using `accept_frac`. The manuscript analyses
generally used a 2.5 percent acceptance fraction before donor balancing.
