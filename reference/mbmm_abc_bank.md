# Fit ABC posteriors for multiple subjects or donors

Applies
[`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md)
to each subject in a named observation list and stacks the accepted
particles.

## Usage

``` r
mbmm_abc_bank(
  obs_list,
  simulate_fn,
  prior_sample_fn,
  n_draws = 10000,
  n_keep = NULL,
  accept_frac = 0.025,
  distance_fn = NULL,
  seed = 1,
  failure = c("discard", "error"),
  verbose = TRUE
)
```

## Arguments

- obs_list:

  Named list of subject observation data.frames.

- simulate_fn:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- prior_sample_fn:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- n_draws:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- n_keep:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- accept_frac:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- distance_fn:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- seed:

  Base random seed. Subject `i` uses `seed + i`.

- failure:

  See
  [`mbmm_abc_fit`](https://linhdo154.github.io/MBMM/reference/mbmm_abc_fit.md).

- verbose:

  Print progress.

## Value

A stacked particle data.frame with donor IDs. Subject-level ABC
diagnostics are stored in the `"abc_diagnostics"` attribute.
