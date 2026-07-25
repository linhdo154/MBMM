# Pool per-donor ABC particles into a mechanistic reference bank

Constructs a balanced mechanistic reference bank by retaining the same
number of accepted ABC particles from every eligible donor.

## Usage

``` r
mbmm_pool_bank(
  particles,
  n_per_donor = NULL,
  donor_col = "id",
  eligible_donors = NULL,
  allow_fewer = FALSE,
  seed = NULL
)
```

## Arguments

- particles:

  A data.frame containing accepted ABC particles, with one row per
  particle.

- n_per_donor:

  Number of particles to retain per donor. If NULL, uses the minimum
  available count across eligible donors.

- donor_col:

  Name of the donor ID column.

- eligible_donors:

  Optional vector of donor IDs to include. This can be used for
  application-specific eligibility rules, such as requiring a minimum
  number of detected measurements.

- allow_fewer:

  If FALSE, stop when a donor has fewer particles than `n_per_donor`. If
  TRUE, retain all available particles for that donor.

- seed:

  Optional random seed.

## Value

A pooled particle bank with a balanced donor contribution when
`allow_fewer = FALSE`.

## Details

When `n_per_donor = NULL`, the smallest available donor-specific
particle count is used, so that every donor contributes equally.

If `n_per_donor` is supplied and some donor has fewer particles than
requested, the default behavior is to stop rather than silently
construct an unbalanced bank. Set `allow_fewer = TRUE` only when unequal
donor-specific particle counts are intentional.
