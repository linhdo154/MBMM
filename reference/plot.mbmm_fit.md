# Plot the bank log-likelihood-ratio diagnostic for an MBMM fit

For a K = 2 MBMM fit, plots each subject's \\\log
L\_{i,\code{classes\[1\]}} - \log L\_{i,\code{classes\[2\]}}\\ (positive
values favor the first class's reference bank), grouped and jittered by
MAP posterior class – the diagnostic used throughout Do, McCarthy &
Rogers to sanity-check that the reference banks separate subjects before
any gate/covariate information is used.

## Usage

``` r
# S3 method for class 'mbmm_fit'
plot(x, ...)
```

## Arguments

- x:

  an `"mbmm_fit"` object from
  [`mbmm`](https://linhdo154.github.io/MBMM/reference/mbmm.md)

- ...:

  additional arguments passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) call

## Value

`x`, invisibly
