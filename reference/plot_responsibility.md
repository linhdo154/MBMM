# Plot posterior responsibilities against simulated truth

Intended for simulations, where the true latent state is known. In real
data there is no truth column to supply: the responsibility is the
inferred probability of the state rather than a calibrated diagnostic
label.

## Usage

``` r
plot_responsibility(fit, truth, title = NULL)
```

## Arguments

- fit:

  An object returned by
  [`fit_gate()`](https://linhdo154.github.io/MBMM/reference/fit_gate.md).

- truth:

  Binary vector of known simulated states (0 or 1).

- title:

  Optional plot title.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
n <- 30; x <- rnorm(n); X <- cbind(1, x)
truth <- rbinom(n, 1, plogis(-0.5 - x))
fit <- fit_gate(rnorm(n, truth), rnorm(n, 1 - truth), X)
plot_responsibility(fit, truth)
```
