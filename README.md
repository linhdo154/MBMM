<!-- README.md is generated from README.Rmd. Please edit that file -->

# mbmm

<!-- badges: start -->
[![R-CMD-check](https://github.com/USERNAME/mbmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USERNAME/mbmm/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/USERNAME/mbmm/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/USERNAME/mbmm/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`mbmm` implements the **Mechanistic Bayesian Mixture Model (MBMM)**
framework of Do, McCarthy, and Rogers (2026), *A Mechanistic Bayesian
Framework for Latent-State Inference in Vaccine and Infection
Immunology*.

MBMM infers posterior probabilities of a latent immunological state
(e.g. productive infection vs. protection/transient detection, or a
boost-responsive vs. weak/stable antibody phenotype) from sparse
longitudinal biomarker trajectories. It does this by embedding an
application-specific ODE trajectory model inside a Bayesian mixture
likelihood evaluated against Approximate Bayesian Computation (ABC)
reference banks, and linking covariates to latent-state probabilities
through a logistic (K = 2) or multinomial (K > 2) gate model.

The package is deliberately domain-agnostic: you supply the
mechanistic trajectory model (ODE simulator), the ABC prior, the
observation model, and the covariate/gate data; `mbmm` supplies the
shared inferential engine — ABC bank construction, marginal
reference-bank likelihoods with leave-one-donor-out scoring,
EM/direct-optimization gate fitting, likelihood-ratio and permutation
testing, and fitted-threshold summaries.

## Installation

You can install the development version of `mbmm` from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("USERNAME/mbmm")
```

## Getting started

The **Introduction to mbmm** vignette walks through the three stages
of a typical analysis — fitting reference banks, scoring subjects, and
fitting the gate model — on simulated correlate-of-protection data. A
second vignette, **Simulation demo**, walks through generating
synthetic cohorts to validate the pipeline end to end. See the
[Articles](articles/index.html) tab for both.

``` r
library(mbmm)
## basic example code
```

## Citation

If you use `mbmm` in your work, please cite:

> Do, T., McCarthy, J., & Rogers, B. (2026). A Mechanistic Bayesian
> Framework for Latent-State Inference in Vaccine and Infection
> Immunology.

## License

MIT © Linh Do
