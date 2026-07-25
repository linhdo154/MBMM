# Safe log10 with a positive floor

Values below `floor`, including zero and negative values, are replaced
by `floor` before transformation. Missing values remain missing.

## Usage

``` r
mbmm_safe_log10(x, floor = 0.01)
```

## Arguments

- x:

  Numeric vector.

- floor:

  Smallest positive value passed to `log10`.

## Value

Numeric vector on the log10 scale.
