# Internal: observed-data marginal log-likelihood for a binary MBMM gate

Computes

## Usage

``` r
.mbmm_obs_ll_binary(beta, X, logL_pos, logL_neg)
```

## Details

\$\$ \ell(\beta) = \sum_i \log\left\[ \pi_i(\beta)L\_{i1} +
\\1-\pi_i(\beta)\\L\_{i0} \right\]. \$\$

This function returns the unpenalized observed-data likelihood.
