# Back-transform an MBMM gate threshold

Converts a threshold from the covariate scale used in the gate model
back to the original biomarker scale.

## Usage

``` r
mbmm_backtransform_threshold(
  threshold,
  center = 0,
  scale = 1,
  transform = c("identity", "log10", "log")
)
```

## Arguments

- threshold:

  Numeric threshold on the fitted gate-covariate scale.

- center:

  Center used when standardizing the covariate.

- scale:

  Scale used when standardizing the covariate.

- transform:

  Transformation applied before standardization. One of `"identity"`,
  `"log10"`, or `"log"`.

## Value

Numeric threshold on the original biomarker scale.

## Details

If the fitted covariate is a standardized variable

\$\$ z = \frac{x - \mu}{s}, \$\$

the model-scale threshold is first converted to

\$\$ x = \mu + sz. \$\$

An optional inverse transformation is then applied. For example, when
the standardized covariate was constructed from `log10(biomarker)`, use
`transform = "log10"`.

## Examples

``` r
mbmm_backtransform_threshold(
  threshold = 1.2,
  center = 2,
  scale = 0.5,
  transform = "log10"
)
#> [1] 398.1072
```
