# mbmm

`mbmm` implements the **Mechanistic Bayesian Mixture Model (MBMM)**
framework of Do, McCarthy, and Rogers (2026), *A Mechanistic Bayesian
Framework for Latent-State Inference in Vaccine and Infection
Immunology*.

MBMM infers posterior probabilities of latent immunological states, such
as productive infection versus non-establishment or transient/abortive
detection, or a boost-responsive versus weak/stable antibody phenotype,
from sparse longitudinal biomarker trajectories.

The framework embeds an application-specific ordinary differential
equation trajectory model inside a Bayesian mixture likelihood evaluated
against Approximate Bayesian Computation (ABC) reference banks.
Covariates are linked to latent-state probabilities through a logistic
gate for two classes or a multinomial gate for more than two classes.

The package is deliberately domain-agnostic. Users supply the
mechanistic trajectory model, ABC prior, observation model, and
covariate data. `mbmm` supplies the shared inferential engine,
including:

- ABC reference-bank construction;
- marginal reference-bank likelihoods;
- leave-one-donor-out scoring;
- logistic and multinomial gate fitting;
- likelihood-ratio and permutation tests;
- fitted-threshold estimation and uncertainty summaries;
- reduced manuscript-style simulation tools.

## Installation

You can install the development version of `mbmm` from GitHub with:

``` r

# install.packages("remotes")
remotes::install_github("linhdo154/MBMM")
```

## Getting started

The [Introduction to
mbmm](https://linhdo154.github.io/MBMM/articles/mbmm-intro.html)
vignette walks through the main stages of an MBMM analysis: constructing
reference banks, scoring subjects, and fitting the covariate-linked gate
model using simulated correlate-of-protection data.

The [reduced manuscript-style
simulation](https://linhdo154.github.io/MBMM/articles/mbmm-simulation-demo.html)
vignette demonstrates how to generate synthetic cohorts and run the
pipeline end to end.

The complete package website is available at:

<https://linhdo154.github.io/MBMM/>

``` r

library(mbmm)
```

## Citation

To obtain the recommended package citation in R, run:

``` r

citation("mbmm")
```

The methodological framework is described in:

> Do, L., McCarthy, G., and Rogers, B. (2026).  
> *A Mechanistic Bayesian Framework for Latent-State Inference in
> Vaccine and Infection Immunology.*

## License

MIT © Linh Do
