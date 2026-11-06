# Screen candidates to the ones participants actually support

A candidate family large enough to cover the observed data is far larger
than the set any one participant is compatible with. Retaining the `nu`
best candidates per participant and pooling them keeps the support the
data speak to and makes the weight estimation tractable.

## Usage

``` r
screen_candidates(loglik, nu = 200L)
```

## Arguments

- loglik:

  A participant-by-candidate log-likelihood matrix.

- nu:

  Number of candidates retained per participant. Larger `nu` keeps more
  of the family at the cost of a bigger weight problem; the retained
  count is not usually binding, because the maximiser on a fixed support
  has at most one positive weight per participant.

## Value

A list with the pooled candidate indices `keep`, the screened matrix
`loglik`, and the number of candidates `before` and `after`.

## Details

The screen uses the measurements, so the pooled support is a function of
the donor data. That is a data-dependent step: it must be run on the
donors only, never on a participant who will later be scored against the
result, or the honesty condition fails. It also means the support is
random, which the fixed-support consistency theory does not cover.

## See also

[`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md),
[`screen_positive()`](https://linhdo154.github.io/MBMM/reference/screen_positive.md)

## Examples

``` r
set.seed(1)
LL <- matrix(rnorm(30 * 200), 30, 200)
s <- screen_candidates(LL, nu = 20)
c(before = s$before, after = s$after)
#> before  after 
#>    200    189 
```
