# Reference-bank class-conditional log-likelihood for one subject

Computes the class-conditional mechanistic likelihood for one subject by
averaging trajectory likelihoods over a class-specific reference bank.

## Usage

``` r
mbmm_class_loglik(
  obs,
  bank_pred,
  sigma,
  exclude_donor = NULL,
  donor_col = "id",
  bundle_col = "bundle_id",
  log_scale = TRUE,
  donor_weighting = c("equal", "particle")
)
```

## Arguments

- obs:

  This subject's observations. A data.frame containing `time` and `y`,
  and optionally `target`, `censored`, and `lod`.

- bank_pred:

  Long-format bank predictions for one class or phenotype. Must contain
  the columns named by `donor_col` and `bundle_col`, together with
  `time`, `pred`, and optionally `target`.

- sigma:

  Measurement-error SD on the Gaussian likelihood scale. Either one
  positive numeric value or, for multitarget data, a named numeric
  vector with one value for every target.

- exclude_donor:

  Optional donor ID or vector of donor IDs to remove before scoring. Use
  this for leave-one-subject/donor-out scoring.

- donor_col:

  Name of the donor column in `bank_pred`.

- bundle_col:

  Name of the within-donor bundle/particle column in `bank_pred`.

- log_scale:

  If `TRUE`, `y`, `pred`, and `lod` are assumed to already be on the
  log10 scale. If `FALSE`, all three are log10-transformed internally.

- donor_weighting:

  Either `"equal"` or `"particle"`. `"equal"` first averages particles
  within donor and then donors, matching the manuscript reference-bank
  construction. `"particle"` averages all donor-particle bundles
  together.

## Value

A scalar class-conditional log-likelihood. Returns `-Inf` when no
eligible bank draw overlaps the subject's observations.

## Details

For equal donor weighting, the function first averages over the retained
particles/bundles contributed by each donor and then averages across
donors:

\$\$ L\_{ik} = \frac{1}{D_k} \sum\_{d=1}^{D_k} \left\[ \frac{1}{M\_{kd}}
\sum\_{m=1}^{M\_{kd}} p(Y_i \mid Z_i = k, \theta\_{kdm}) \right\]. \$\$

This two-stage averaging ensures that a donor with more retained
particles does not receive greater weight in the reference-bank
likelihood.

Detected observations contribute a Gaussian density. Left-censored
observations contribute the Gaussian probability of falling below the
LOD.

## Examples

``` r
obs <- data.frame(
  time = c(0, 1, 2),
  y = c(1, 2, 3)
)

close_bank <- data.frame(
  id = "donorA",
  bundle_id = 1,
  time = c(0, 1, 2),
  pred = c(1.05, 1.95, 3.1)
)

far_bank <- data.frame(
  id = "donorB",
  bundle_id = 1,
  time = c(0, 1, 2),
  pred = c(5, 5, 5)
)

mbmm_class_loglik(obs, close_bank, sigma = 0.2)
#> [1] 1.883998
mbmm_class_loglik(obs, far_bank, sigma = 0.2)
#> [1] -360.4285
```
