# Join target-specific particle banks into complete multitarget bundles

Constructs joint mechanistic bundles for multitarget applications. For
each donor and each bundle, one accepted particle is sampled
independently from every target-specific donor pool.

## Usage

``` r
mbmm_join_targets(
  particles_by_target,
  n_bundles = 80,
  donor_col = "id",
  replace = TRUE,
  seed = NULL
)
```

## Arguments

- particles_by_target:

  Named list of particle data.frames, one per target.

- n_bundles:

  Number of joint bundles to construct per donor.

- donor_col:

  Name of the donor ID column.

- replace:

  Whether to sample particles with replacement within each
  target-specific donor pool. The SDY2354 construction typically uses
  replacement because the requested number of joint bundles may exceed
  the number of retained particles for an individual target.

- seed:

  Optional random seed.

## Value

A long data.frame with one row per donor, bundle, and target.

## Details

Only donors represented in every target are eligible. Every returned
donor-bundle combination is guaranteed to contain exactly one row for
every target.
