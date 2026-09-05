# Changelog

## metaDyn 1.0.4

### Patch

- Improved convergence assessment and Hessian-rescue behavior, including
  more robust selection of the best candidate fit.
- Treats lower-bound solutions for softplus covariance-diagonal
  parameters as admissible rather than automatically classifying them as
  convergence failures.
- Hardened parameter-bound detection, nudging, and bound relaxation.

## metaDyn 1.0.3

CRAN release: 2026-08-02

### Patch

- Revised the model-implied mean vector and covariance matrix for
  distal-outcome models so that distal outcomes are predicted by latent
  effect sizes rather than observed effect-size estimates.
- Added support for treating covariates as either fixed definition
  variables or stochastic variables modeled jointly with effect sizes
  and distal outcomes.
- Added Monte Carlo method confidence intervals.

## metaDyn 1.0.1

CRAN release: 2026-03-24

### Patch

- Initial CRAN release.

## metaDyn 1.0.0

### Patch

- Initial stable version.
