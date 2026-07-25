# Internal: observed-data marginal log-likelihood for a multinomial gate

The final latent class is the reference class, with its coefficient
vector fixed at zero.

## Usage

``` r
.mbmm_obs_ll_multinom(beta_vec, X, logL_mat, K)
```

## Arguments

- beta_vec:

  Numeric vector containing coefficients for the first `K - 1` classes.

- X:

  Gate design matrix.

- logL_mat:

  Subject-by-class matrix of class log-likelihoods.

- K:

  Number of latent classes.
