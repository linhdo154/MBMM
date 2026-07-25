# Internal: ridge penalty for a multinomial gate

Applies an L2 penalty to all non-intercept coefficients for the first
`K - 1` classes.

## Usage

``` r
.mbmm_ridge_penalty_multinom(beta_vec, p, K, ridge)
```
