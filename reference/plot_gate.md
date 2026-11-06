# Plot a fitted covariate gate

Displays the fitted prior state probability against one covariate. For a
simulation, optional true states and generating probabilities make
recovery visible. For real data, omit both and the fitted gate alone is
shown.

## Usage

``` r
plot_gate(
  fit,
  x,
  truth = NULL,
  true_probability = NULL,
  xlab = "Covariate",
  title = NULL
)
```

## Arguments

- fit:

  An object returned by
  [`fit_gate()`](https://linhdo154.github.io/MBMM/reference/fit_gate.md).

- x:

  Numeric covariate used in the gate.

- truth:

  Optional binary vector of known simulated latent states.

- true_probability:

  Optional vector of generating state probabilities.

- xlab:

  Horizontal-axis label.

- title:

  Optional plot title.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
n <- 40; x <- rnorm(n); prob <- plogis(-0.5 - x)
truth <- rbinom(n, 1, prob); X <- cbind(1, x)
fit <- fit_gate(rnorm(n, truth), rnorm(n, 1 - truth), X)
plot_gate(fit, x, truth, prob)
```
