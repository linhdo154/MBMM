# Restrict a bank to a region of parameter space (Route II)

Defines a latent class as a prespecified region of kinetic-parameter
space and recovers its kinetic distribution by restricting an estimated
marginal population to that region and renormalising. Because the class
prevalence cancels, the population used to build the bank need not
consist only of class members, nor share the analysis population's class
prevalence.

## Usage

``` r
restrict_bank(bank, in_region)
```

## Arguments

- bank:

  An `mbmm_bank` from
  [`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md),
  or a numeric weight vector.

- in_region:

  Logical vector, one entry per candidate, `TRUE` for candidates inside
  the region.

## Value

A list with the renormalised `weights`, the `mass` the region carried in
the marginal, and `in_region`.

## Details

The region must be fixed before any outcome is observed, and must carry
population mass bounded away from zero: restriction divides by that
mass, so a thin region amplifies any error in the marginal.

## References

Result D of the accompanying manuscript.

## Examples

``` r
set.seed(1)
b <- fit_bank(matrix(rnorm(20 * 6), 20, 6))
restrict_bank(b, in_region = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))$mass
#> [1] 0.4237955
```
