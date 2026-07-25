# Fit the MBMM gate model

Fits the covariate-linked gate model using fixed class-conditional
mechanistic log-likelihoods.

## Usage

``` r
mbmm_fit_gate(
  logL,
  formula,
  data,
  classes = NULL,
  ridge = 1e-04,
  penalty = "ridge",
  n_perm = 0,
  seed = NULL,
  id_col = NULL,
  maxit = 4000
)
```

## Arguments

- logL:

  Data.frame or matrix of class-conditional log-likelihoods. Columns are
  normally named `logL_<class>`.

- formula:

  One-sided formula defining the gate covariates, such as
  `~ baseline_antibody` or `~ arm`.

- data:

  Data.frame containing gate covariates.

- classes:

  Character vector defining the class order. For two classes,
  `classes[1]` is the positive/event class and `classes[2]` is the
  reference class. For more than two classes, the final class is the
  multinomial reference class.

- ridge:

  Nonnegative L2 penalty applied to non-intercept coefficients when
  `penalty = "ridge"`.

- penalty:

  Either `"ridge"` or `"none"`.

- n_perm:

  Number of permutations for an optional permutation p-value. The
  implemented procedure permutes the complete full-model design matrix
  and is appropriate only when the tested non-intercept covariates are
  exchangeable as a block.

- seed:

  Optional random seed for permutation testing.

- id_col:

  Optional subject-ID column in `data`, copied to the posterior table.

- maxit:

  Maximum number of optimizer iterations.

## Value

An object of class `"mbmm_gate"`.

## Details

For a binary latent state, the observed-data marginal likelihood is

\$\$ \ell(\beta) = \sum_i \log\left\[ \pi_i(\beta)L\_{i1} +
\\1-\pi_i(\beta)\\L\_{i0} \right\]. \$\$

The parameters are estimated by direct optimization of this marginal
likelihood. An optional ridge penalty can be applied to non-intercept
coefficients for numerical stabilization. Ridge is not part of the core
MBMM model.

The reported likelihood-ratio statistic is evaluated using the
unpenalized observed-data marginal likelihood at the fitted full- and
null-model coefficient estimates:

\$\$ 2\\\ell(\widehat{\beta}\_{full}) - \ell(\widehat{\beta}\_{null})\\.
\$\$

## Examples

``` r
set.seed(1)

n <- 150
x <- rnorm(n)

z <- rbinom(
  n,
  size = 1,
  prob = mbmm_inv_logit(-0.5 + 1.5 * x)
)

logL <- data.frame(
  logL_pos = ifelse(
    z == 1,
    log(0.98),
    log(0.02)
  ),
  logL_neg = ifelse(
    z == 1,
    log(0.02),
    log(0.98)
  )
)

fit <- mbmm_fit_gate(
  logL = logL,
  formula = ~x,
  data = data.frame(x = x),
  classes = c("pos", "neg"),
  penalty = "ridge",
  ridge = 1e-4
)

fit
#> MBMM gate model
#>   classes: pos vs neg 
#>   method: direct-marginal-optimization 
#>   penalty: ridge 
#>   ridge: 1e-04 
#>   n subjects: 150 
#>   optimizer convergence: full = 0 , null = 0 
#>   unpenalized logLik: full = -89.1592 , null = -102.353 
#>   LRT = 26.39 on df = 1  -> p = 2.794e-07 
#> 
#> Gate coefficients:
#>         term   estimate odds_ratio
#>  (Intercept) -0.3977403  0.6718365
#>            x  1.1465444  3.1472984
```
