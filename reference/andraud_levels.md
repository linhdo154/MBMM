# Class coordinate for a prime-boost question

Returns the log-ratio that defines a class region. With `ref = "C1"` the
coordinate is \\g = \log\_{10}(D_1/C_1)\\, asking whether the *boost*
added a durable gain; with `ref = "A0"` it is \\g =
\log\_{10}(D_1/A_0)\\, asking whether the *regimen* left one. In both
cases class 1 is \\\\g \> 0\\\\, so zero is the only boundary and no
threshold is imposed.

## Usage

``` r
andraud_levels(A0, C1, D1, ref = c("C1", "A0"))
```

## Arguments

- A0, C1, D1:

  Baseline, prime-only plateau, post-boost plateau.

- ref:

  Either `"C1"` or `"A0"`.

## Value

A numeric vector of \\g\\.

## Details

The reference level belongs to the candidate design and must be fixed
before fitting. Two questions therefore need two fitted populations;
they are not two regions carved from one.

## Examples

``` r
andraud_levels(A0 = 400, C1 = 1200, D1 = 3000, ref = "C1")
#> [1] 0.39794
andraud_levels(A0 = 400, C1 = 1200, D1 = 3000, ref = "A0")
#> [1] 0.8750613
```
