# mbmm

<!-- badges: start -->
[![R-CMD-check](https://github.com/linhdo154/MBMM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/linhdo154/MBMM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/linhdo154/MBMM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/linhdo154/MBMM/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

**mbmm** implements the **Mechanistic Bayes Mixture Model (MBMM)**
framework described by Do, McCarthy, and Rogers (2026) in *Mechanistic
Bayes Mixture Models*.

MBMM estimates how a covariate changes the probability of a latent biological
state expressed through a longitudinal mechanistic trajectory. Examples
include productive infection versus non-establishment, or a positive durable
antibody gain versus no durable gain.

Rather than reducing a trajectory to a thresholded label, MBMM evaluates each
participant against populations of mechanistic trajectories. A
covariate-dependent logistic gate supplies the prior state probability, and
the trajectory evidence and gate combine into a posterior responsibility.
Participants whose trajectories do not settle the state contribute according
to that uncertainty instead of being forced into a class.

The package is model-agnostic. Users may supply any positive trajectory
function with arguments **time** and **theta**, including a custom ODE solved
with **deSolve**. Built-in eclipse-phase viral and prime-boost antibody models
are included for the manuscript examples.

The current package provides:

- censoring-aware participant-by-trajectory likelihoods;
- integration over latent timing offsets and other nuisance coordinates;
- entropy-penalised kinetic-population estimation;
- independent-donor and exact leave-one-out bank construction;
- Route I state-specific reference banks;
- Route II restriction to a prespecified parameter-space region;
- logistic gate fitting, including optional Firth penalisation;
- likelihood-ratio or penalised-objective statistics with permutation tests;
- diagnostic plots for banks, responsibilities, and fitted gates; and
- fully synthetic examples with known truth and executable recovery checks.

## Installation

Install the development version from GitHub:

~~~r
# install.packages("remotes")
remotes::install_github("linhdo154/MBMM")
~~~

Then load the package:

~~~r
library(mbmm)
~~~

## Getting started

The package has four linked vignettes:

- [Should I use MBMM?](https://linhdo154.github.io/MBMM/articles/mbmm.html)
  explains suitable scientific questions, required assumptions, core
  functions, diagnostics, and reporting conditions.
- [Bring your own ODE](https://linhdo154.github.io/MBMM/articles/custom-ode.html)
  reconstructs the built-in models from scratch, shows both familiar
  parameter boxes and identifiable bank designs, and fits a custom
  constant-beta HIV model.
- [Route I](https://linhdo154.github.io/MBMM/articles/pipeline.html)
  uses an independent pure reference bank in a productive-infection
  simulation and checks recovery against known truth.
- [Route II](https://linhdo154.github.io/MBMM/articles/sdy2354.html)
  defines a durable-gain state as a prespecified parameter region, uses exact
  leave-one-out banks, and checks recovery against known truth.

The complete package website is available at:

<https://linhdo154.github.io/MBMM/>

## Minimal workflow

~~~r
# Score each participant against each candidate trajectory
LL <- trajectory_loglik(
  obs, candidates, trajectory = my_model, sigma = 0.3
)

# Estimate a kinetic population from an appropriate donor likelihood matrix
bank <- fit_bank(LL_donors, lambda = 1)

# Form the two class likelihoods
l1 <- class_loglik(LL, bank$weights)
l0 <- background_loglik(obs, sigma = 0.3)

# Fit and test the covariate gate
fit <- fit_gate(l1, l0, X)
test <- gate_test(l1, l0, X, B = 999)
~~~

The bank used to score a participant must not have been fitted or screened
using that participant's measurements. Use a separate donor sample or exact
cross-fitting with **bank_loo()**.

## Data and reproduction

No participant-level or packaged SDY2354 dataset is distributed with mbmm.
To reproduce the manuscript's SDY2354 analysis, obtain SDY2354 directly from
ImmPort under its user agreement and provide the prepared table to
**inst/reproduce/sdy2354_table3.R**.

**simulate_mbmm()** supplies a fully synthetic Route I cohort with known
latent states.

## Citation

To obtain the recommended package citation in R:

~~~r
citation("mbmm")
~~~

The methodological framework is described in:

> Do, L., McCarthy, J., and Rogers, B. (2026).  
> *Mechanistic Bayes Mixture Models.*

## License

MIT © Linh Do
