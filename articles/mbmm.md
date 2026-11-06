# Should I use MBMM? A scientific and practical guide

``` r

library(mbmm)
```

## Start with the scientific question

MBMM is designed for a question of the form:

> **Does a covariate change the probability of a latent biological state
> that is expressed through a longitudinal mechanistic trajectory?**

The state is not a thresholded measurement. It is a biological condition
such as productive infection or a durable antibody gain. The
longitudinal record contains evidence about that state, and a
mechanistic model says what records from each state should look like.

MBMM separates two tasks:

1.  the **trajectory likelihood** asks what state is compatible with the
    participant’s longitudinal measurements; and
2.  the **gate** asks how a covariate changes the probability of that
    state.

The primary scientific estimand is the gate coefficient. The trajectory
bank is the mechanism used to form state evidence; it is not usually the
endpoint of the study.

## When MBMM is a good match

MBMM is potentially useful when all of the following are true:

- the endpoint is a biologically meaningful but unobserved state;
- state membership changes the shape, scale, timing, or durability of a
  longitudinal trajectory;
- a mechanistic model can generate plausible trajectories for that
  process;
- repeated measurements contain more information than a single
  thresholded label;
- the covariate belongs in the state-probability model, not in the
  trajectory likelihood used to infer the state;
- a state-specific kinetic population can be constructed by Route I or
  Route II; and
- the bank used to score a participant can be built without that
  participant’s measurements.

Examples include:

- whether baseline antibody changes the probability of productive
  infection, inferred from a complete viral-load trajectory;
- whether vaccine regimen changes the probability of a durable
  post-boost gain, inferred from prime-boost antibody kinetics; and
- whether a treatment changes the probability of entering a
  mechanistically defined response phenotype rather than merely crossing
  a visit-specific threshold.

## When MBMM is not the right model

Do not use MBMM merely because the data are longitudinal. A different
model is usually preferable when:

- the endpoint is directly observed without meaningful classification
  error;
- the scientific target is a continuous kinetic parameter rather than a
  latent class;
- there is no defensible mechanistic model or no longitudinal
  information capable of separating the states;
- the class boundary was chosen after looking at the same outcomes and
  cannot be cross-fitted;
- neither a credible state-specific reference population nor a
  prespecified identifiable parameter-space region is available;
- the covariate changes kinetics *within* a state, violating the
  gate-only assumption, unless that dependence is explicitly modelled;
  or
- the donor kinetic population is not transportable to the analysis
  cohort.

An honest computation cannot repair a scientifically invalid state
definition, a non-transportable reference population, or an inadequate
mechanistic model.

## Translate the question into MBMM objects

Before writing code, fill in this table for the proposed analysis.

| Object | Question to answer |
|----|----|
| Longitudinal data $`Y_i`$ | Which repeated measurements carry evidence about the state? |
| Latent state $`Z_i`$ | What biological condition does state 1 mean? What does state 0 mean? |
| Mechanistic trajectory $`\mu(t; \theta)`$ | Which kinetic parameters generate the observed shape? |
| Covariates $`x_i`$ | Which variables are hypothesised to shift state probability? |
| Gate $`pi_i`$ | What probability is being modelled before using the trajectory evidence? |
| Kinetic population $`G_k`$ | What population of trajectories represents each state? |
| Class likelihood $`L_{ik}`$ | How will trajectory evidence be integrated over that population? |
| Responsibility $`w_i`$ | What participant-level posterior state probability will be reported? |

If any row has no scientifically interpretable answer, the analysis is
not ready to fit.

## Choose Route I or Route II

### Route I: known state-specific reference population

Use Route I when a separate donor sample can credibly represent state 1.

Required conditions:

- **Purity:** donors represent the intended state, not an unrestricted
  mixture of states.
- **Transportability:** the within-state kinetic distribution applies to
  the analysis cohort.
- **Honesty:** the complete donor-bank construction is independent of
  the analysis participants.
- **Coverage:** the candidate family covers the observed analysis
  trajectories without relying on implausible boundary candidates.

Continue to the [Route I
simulation](https://linhdo154.github.io/MBMM/articles/pipeline.md),
which builds an independent productive-infection bank with the
eclipse-phase viral ODE and checks recovery against known truth.

``` r

vignette("pipeline", package = "mbmm")
```

### Route II: prespecified region of kinetic parameter space

Use Route II when no labelled donor sample exists but the scientific
state is a region $`S_1`$ fixed before fitting, such as a positive
durable gain.

Required conditions:

- **Target validity:** the region corresponds to the scientific state
  and the boundary is not selected to optimise the observed result.
- **Identifiability:** the longitudinal measurements distinguish the
  relevant kinetic coordinate.
- **Region support:** both $`S_1`$ and its complement carry adequate
  fitted marginal mass; restriction by nearly zero mass is unstable.
- **Honesty:** every outcome-adaptive step, including support screening
  and weight fitting, excludes the participant being scored.
- **Exchangeability:** the remaining participants are a suitable
  population from which to estimate the held-out participant’s marginal
  kinetics.

Continue to the [Route II
simulation](https://linhdo154.github.io/MBMM/articles/sdy2354.md), which
defines durable gain by $`g=log_{10}(D_1/C_1)>0`$, fits exact
leave-one-out marginal banks, and checks recovery against known truth.

``` r

vignette("sdy2354", package = "mbmm")
```

## Core package functions, in scientific order

MBMM is not limited to its built-in models. Any positive trajectory
function with arguments **time** and **theta** can be passed to
**trajectory_loglik()**. The [custom ODE
walkthrough](https://linhdo154.github.io/MBMM/articles/custom-ode.md)
reconstructs both built-in models from scratch, explains raw parameter
boxes and identifiable bank designs, and fits a new HIV constant-beta
model.

| Scientific step | MBMM functions | What to verify |
|----|----|----|
| Generate mechanistic trajectories | [`eclipse_trajectory()`](https://linhdo154.github.io/MBMM/reference/eclipse_trajectory.md), [`andraud_trajectory()`](https://linhdo154.github.io/MBMM/reference/andraud_trajectory.md) | Curves are biologically plausible over a dense time grid |
| Construct a candidate family | [`eclipse_design()`](https://linhdo154.github.io/MBMM/reference/eclipse_design.md), user-defined candidate data frames | Bounds and screens are fixed independently of scored outcomes |
| Exclude invalid trajectories | [`screen_positive()`](https://linhdo154.github.io/MBMM/reference/screen_positive.md) | Screening rule is mechanistic, or is repeated inside cross-fitting if data-adaptive |
| Score observations against candidates | [`trajectory_loglik()`](https://linhdo154.github.io/MBMM/reference/trajectory_loglik.md), [`obs_loglik()`](https://linhdo154.github.io/MBMM/reference/obs_loglik.md) | Noise, censoring, visit times, and any timing uncertainty match the data-generating assumptions |
| Integrate nuisance coordinates | [`pool_marginal()`](https://linhdo154.github.io/MBMM/reference/pool_marginal.md), [`jitter_marginal()`](https://linhdo154.github.io/MBMM/reference/jitter_marginal.md) | Nuisance variation is integrated rather than allowed to redefine the class |
| Estimate a kinetic population | [`fit_bank()`](https://linhdo154.github.io/MBMM/reference/fit_bank.md) | Donors are pure/appropriate and independent of scored participants |
| Cross-fit an in-sample population | [`bank_loo()`](https://linhdo154.github.io/MBMM/reference/bank_loo.md) | All adaptive steps exclude the held-out participant |
| Define a Route II class | [`restrict_bank()`](https://linhdo154.github.io/MBMM/reference/restrict_bank.md) | Region is prespecified and has adequate mass |
| Form a class likelihood | [`class_loglik()`](https://linhdo154.github.io/MBMM/reference/class_loglik.md), [`background_loglik()`](https://linhdo154.github.io/MBMM/reference/background_loglik.md) | $`L_{i1}`$ and $`L_{i0}`$ represent the intended states |
| Fit the covariate gate | [`fit_gate()`](https://linhdo154.github.io/MBMM/reference/fit_gate.md) | Design matrix matches the scientific contrast; separation is handled appropriately |
| Test a gate effect | [`gate_test()`](https://linhdo154.github.io/MBMM/reference/gate_test.md) | Permutation scheme matches exchangeability; penalised statistics are not called likelihood ratios |
| Visualise diagnostics | [`plot_bank()`](https://linhdo154.github.io/MBMM/reference/plot_bank.md), [`plot_responsibility()`](https://linhdo154.github.io/MBMM/reference/plot_responsibility.md), [`plot_gate()`](https://linhdo154.github.io/MBMM/reference/plot_gate.md) | Plots are interpreted as diagnostics, not proof of assumptions |

## Conditions to check before trusting the output

### 1. State meaning and class construction

Write the biological definition of both states without referring to the
fitted result. Route I needs a state-specific reference population.
Route II needs a boundary fixed before the outcome-adaptive fit.

### 2. Observation model

Confirm the measurement scale, residual noise, censoring mechanism,
visit times, and any latent timing offset. The residual scale controls
how sharply likelihood evidence separates the states; choosing it too
small can harden responsibilities artificially.

### 3. Candidate coverage and positivity

Plot candidate trajectories on a dense grid. Check that observed records
are covered, biologically invalid curves are excluded, and fitted weight
does not pile up at arbitrary design boundaries.

### 4. Honesty

The participant being scored must not contribute to the bank, support
screen, noise tuning, or any other outcome-adaptive step used for that
score. Use a fully independent bank or exact cross-fitting.

### 5. Purity or region validity

Route I requires the bank to represent the intended state. Route II
requires the marginal population to be restricted to an identifiable,
prespecified region. Honesty and purity are separate conditions:
leave-one-out fitting can be honest while still estimating the wrong
mixture.

### 6. Transportability or exchangeability

For Route I, ask whether donor kinetics within the state apply to the
analysis cohort. For Route II, ask whether the other participants are
suitable donors for the held-out participant’s marginal kinetic
distribution.

### 7. Stability and region mass

Report effective support, sensitivity to candidate design and residual
scale, and, for Route II, the mass of every restricted region in every
fold. These are stability diagnostics; because the oracle bank is
unknown, they are not verified bounds on oracle error.

### 8. Gate interpretation and inference

The fitted gate probability `fit$pi` uses covariates before the
participant’s trajectory. The responsibility `fit$w` combines that prior
probability with the class likelihood ratio. The coefficient `fit$gamma`
is a log-odds association with latent-state probability, conditional on
the model.

Use permutation inference when the reference distribution of the plug-in
statistic is not established. With `firth = TRUE`,
[`gate_test()`](https://linhdo154.github.io/MBMM/reference/gate_test.md)
reports a penalised-objective difference and a permutation p-value; it
is not an ordinary likelihood-ratio test.

## A minimum reporting checklist

Report at least:

- the biological definitions of states 0 and 1;
- Route I or Route II and why that route is valid;
- the mechanistic and observation models;
- candidate bounds, screens, integration grids, and residual scale;
- the complete honesty/cross-fitting construction;
- donor purity and transportability, or Route II region definition and
  mass;
- effective bank support and sensitivity analyses;
- responsibility distributions and the fitted gate contrast;
- the exact permutation procedure and number of permutations; and
- limitations, especially within-class covariate effects and kinetic
  mismatch.

## Choose the detailed walkthrough

- [Bring your own ODE: custom mechanistic
  models](https://linhdo154.github.io/MBMM/articles/custom-ode.md)
- [Route I: independent reference-bank
  simulation](https://linhdo154.github.io/MBMM/articles/pipeline.md)
- [Route II: prespecified-region
  simulation](https://linhdo154.github.io/MBMM/articles/sdy2354.md)

Both simulations expose the latent truth, show how each line of code
maps to the MBMM hierarchy, and contain executable recovery checks.
Passing those examples verifies the implementation in those controlled
settings; it does not substitute for checking the scientific conditions
above in a new study.
