# mbmm: Mechanistic Bayesian Mixture Models for Latent Kinetic States

Estimates how a covariate shifts the probability of a latent biological
state that is expressed through a longitudinal trajectory rather than
observed directly. Each state is represented by a population of
mechanistic trajectories (a "bank"); integrating a censoring-aware
trajectory likelihood over that population gives a class likelihood, and
a covariate-dependent gate supplies the prior state probability. An EM
algorithm estimates the gate by weighted logistic regression on the
resulting posterior responsibilities, so participants whose state is
uncertain contribute at the weight the evidence warrants. Includes the
eclipse-phase target-cell-limited viral model and the reduced Andraud
prime-boost antibody model used in the accompanying manuscript.

## See also

Useful links:

- <https://github.com/linhdo154/MBMM>

- <https://linhdo154.github.io/MBMM/>

- Report bugs at <https://github.com/linhdo154/MBMM/issues>

## Author

**Maintainer**: Linh Do <thuylinh.do@duke.edu>

Authors:

- Janice McCarthy

- Bruce Rogers
