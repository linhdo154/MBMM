# Plot a fitted MBMM gate curve against one continuous covariate

For a binary (K = 2) gate fit with a single continuous covariate (plus
optional intercept), draws the fitted gate probability curve \\P(Z =
\code{classes\[1\]} \mid x)\\ over the observed range of that covariate,
with subjects' posterior class probabilities overlaid as points
(jittered on the y-axis at their MAP class for readability). This is the
plot behind Fig. 4/Fig. 6-style panels in Do, McCarthy & Rogers.

## Usage

``` r
# S3 method for class 'mbmm_gate'
plot(x, focal = NULL, hold_at = list(), n_grid = 200, ...)
```

## Arguments

- x:

  an `"mbmm_gate"` object from
  [`mbmm_fit_gate`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md)

- focal:

  name of the continuous covariate to plot on the x-axis. If NULL
  (default), the first non-intercept term in the gate formula is used.

- hold_at:

  named list of fixed values for any other covariates in the gate
  formula (only needed if there is more than one covariate)

- n_grid:

  number of points used to draw the fitted curve

- ...:

  additional arguments passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) call (e.g.
  `main`, `col`)

## Value

the fitted gate object, invisibly
