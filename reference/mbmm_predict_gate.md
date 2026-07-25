# Predict fitted MBMM gate probabilities

Predicts the gate probabilities implied by the fitted covariate model
before incorporating subject-specific longitudinal mechanistic
likelihoods. These are therefore prior latent-class probabilities under
the MBMM gate.

## Usage

``` r
mbmm_predict_gate(fit, newdata)
```

## Arguments

- fit:

  An object of class `"mbmm_gate"` returned by
  [`mbmm_fit_gate`](https://linhdo154.github.io/MBMM/reference/mbmm_fit_gate.md).

- newdata:

  A data.frame containing the variables used in `fit$formula`.

## Value

For a binary gate, a numeric vector giving the probability of
`fit$classes[1]`. For a multinomial gate, a matrix with one column per
latent class.
