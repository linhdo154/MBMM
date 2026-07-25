# Score every subject against class-specific reference banks

Builds the subject-by-class matrix of class-conditional log-likelihoods
used by the MBMM gate model.

## Usage

``` r
mbmm_score_subjects(
  obs_list,
  bank_pred_list,
  sigma,
  donor_col = "id",
  bundle_col = "bundle_id",
  log_scale = TRUE,
  donor_weighting = c("equal", "particle"),
  loo = TRUE,
  verbose = TRUE
)
```

## Arguments

- obs_list:

  Named list of subject observation data.frames.

- bank_pred_list:

  Named list of reference-bank prediction tables, one per bank-based
  class.

- sigma:

  Measurement-error SD passed to
  [`mbmm_class_loglik`](https://linhdo154.github.io/MBMM/reference/mbmm_class_loglik.md).

- donor_col:

  Name of the donor column in each bank table.

- bundle_col:

  Name of the bundle column in each bank table.

- log_scale:

  Whether observations, predictions, and LODs are already on the log10
  likelihood scale.

- donor_weighting:

  Either `"equal"` or `"particle"`.

- loo:

  Logical scalar or named logical vector. If a scalar, the same LOSO
  rule is applied to every class. If named, names must correspond to
  `bank_pred_list`.

- verbose:

  Print progress every ten subjects.

## Value

A data.frame with one row per subject and one `logL_<class>` column per
bank.

## Details

Leave-one-subject/donor-out scoring is controlled by `loo`; it is not
applied automatically unless requested. This allows the same function to
reproduce application-specific choices:

- SDY2354: typically `loo = TRUE`;

- COVID5001 productive bank: use the setting specified by the final
  manuscript analysis.
