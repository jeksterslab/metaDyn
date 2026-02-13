# Fit Multivariate Meta-Analysis (metaVAR via OpenMx)

This function estimates fixed-, random-, or mixed-effects meta-analytic
parameters using per-individual coefficient estimates and their sampling
variance-covariance matrices. Optionally, it fits distal-outcome models
in which between-person outcomes are regressed on between-person
covariates and the meta-analyzed parameters/effect sizes. This function
uses the estimated coefficients and sampling variance-covariance matrix
from each individual fitted using the
[`fitVARMxID::FitVARMxID()`](https://github.com/jeksterslab/fitVARMxID/reference/FitVARMxID.html)
function.

## Usage

``` r
MetaVARMx(
  object,
  x = NULL,
  z = NULL,
  random = TRUE,
  alpha_free = NULL,
  alpha_values = NULL,
  alpha_lbound = NULL,
  alpha_ubound = NULL,
  tau_sqr_diag = FALSE,
  tau_sqr_d_free = NULL,
  tau_sqr_d_values = NULL,
  tau_sqr_d_lbound = NULL,
  tau_sqr_d_ubound = NULL,
  tau_sqr_l_free = NULL,
  tau_sqr_l_values = NULL,
  tau_sqr_l_lbound = NULL,
  tau_sqr_l_ubound = NULL,
  i_sqr_univariate = FALSE,
  gamma_free = NULL,
  gamma_values = NULL,
  gamma_lbound = NULL,
  gamma_ubound = NULL,
  kappa_free = NULL,
  kappa_values = NULL,
  kappa_lbound = NULL,
  kappa_ubound = NULL,
  phi_values = NULL,
  phi_free = NULL,
  phi_lbound = NULL,
  phi_ubound = NULL,
  omega_values = NULL,
  omega_free = NULL,
  omega_lbound = NULL,
  omega_ubound = NULL,
  psi_diag = TRUE,
  psi_d_free = NULL,
  psi_d_values = NULL,
  psi_d_lbound = NULL,
  psi_d_ubound = NULL,
  psi_l_free = NULL,
  psi_l_values = NULL,
  psi_l_lbound = NULL,
  psi_l_ubound = NULL,
  check_estimates = TRUE,
  effects = TRUE,
  set_point = FALSE,
  int_meas = FALSE,
  int_dyn = FALSE,
  cov_meas = FALSE,
  cov_dyn = FALSE,
  diag_cov = "var",
  converged = TRUE,
  vanishing_theta = TRUE,
  theta_tol = 0.001,
  robust_v = FALSE,
  robust = FALSE,
  lb = FALSE,
  alpha = 0.05,
  tries_explore = 100,
  tries_local = 100,
  max_attempts = 10,
  grad_tol = 0.01,
  hess_tol = 1e-08,
  eps = 1e-06,
  factor = 10,
  abs_bnd_tol = 1e-06,
  rel_bnd_tol = 1e-04,
  silent = FALSE,
  seed = NULL,
  ncores = NULL
)
```

## Arguments

- object:

  Output of the
  [`fitVARMxID::FitVARMxID()`](https://github.com/jeksterslab/fitVARMxID/reference/FitVARMxID.html)
  function.

- x:

  An optional list. Each element of the list is a numeric vector of
  covariates.

- z:

  An optional list. Each element of the list is a numeric vector of
  distal outcomes.

- random:

  Logical. If `random = TRUE`, estimates random effects. If
  `random = FALSE`, `tau_sqr` is a null matrix.

- alpha_free:

  Logical vector. Optional vector of free (`TRUE`) parameters for
  `alpha`.

- alpha_values:

  Numeric vector. Optional vector of starting values for `alpha`.

- alpha_lbound:

  Numeric vector. Optional vector of lower bound values for `alpha`.

- alpha_ubound:

  Numeric vector. Optional vector of upper bound values for `alpha`.

- tau_sqr_diag:

  Logical. If `tau_sqr_diag = TRUE`, `tau_sqr` is a diagonal matrix. If
  `tau_sqr_diag = FALSE`, `tau_sqr` is a symmetric matrix.

- tau_sqr_d_free:

  Logical vector indicating free/fixed status of the elements of
  `tau_sqr_d`. If `NULL`, all element of `tau_sqr_d` are free.

- tau_sqr_d_values:

  Numeric vector with starting values for `tau_sqr_d`. If `NULL`,
  defaults to a vector of ones.

- tau_sqr_d_lbound:

  Numeric vector with lower bounds for `tau_sqr_d`. If `NULL`, no lower
  bounds are set.

- tau_sqr_d_ubound:

  Numeric vector with upper bounds for `tau_sqr_d`. If `NULL`, no upper
  bounds are set.

- tau_sqr_l_free:

  Logical matrix indicating which strictly-lower-triangular elements of
  `tau_sqr_l` are free. Ignored if `tau_sqr_diag = TRUE`.

- tau_sqr_l_values:

  Numeric matrix of starting values for the strictly-lower-triangular
  elements of `tau_sqr_l`. If `NULL`, defaults to a null matrix.

- tau_sqr_l_lbound:

  Numeric matrix with lower bounds for `tau_sqr_l`. If `NULL`, no lower
  bounds are set.

- tau_sqr_l_ubound:

  Numeric matrix with upper bounds for `tau_sqr_l`. If `NULL`, no upper
  bounds are set.

- i_sqr_univariate:

  Logical. If `i_sqr_univariate = TRUE`, use the univariate formula for
  \\I^2\\. If `i_sqr_univariate = FALSE`, use the multivariate formula
  for \\I^2\\.

- gamma_free:

  Logical matrix. Optional matrix of free (`TRUE`) parameters for
  `gamma`.

- gamma_values:

  Numeric matrix. Optional matrix of starting values for `gamma`.

- gamma_lbound:

  Numeric matrix. Optional matrix of lower bound values for `gamma`.

- gamma_ubound:

  Numeric matrix. Optional matrix of upper bound values for `gamma`.

- kappa_free:

  Logical vector. Optional vector of free (`TRUE`) parameters for
  `kappa`.

- kappa_values:

  Numeric vector. Optional vector of starting values for `kappa`.

- kappa_lbound:

  Numeric vector. Optional vector of lower bound values for `kappa`.

- kappa_ubound:

  Numeric vector. Optional vector of upper bound values for `kappa`.

- phi_values:

  Numeric matrix. Optional matrix of starting values for `phi`.

- phi_free:

  Logical matrix. Optional matrix of free (`TRUE`) parameters for `phi`.

- phi_lbound:

  Numeric matrix. Optional matrix of lower bound values for `phi`.

- phi_ubound:

  Numeric matrix. Optional matrix of upper bound values for `phi`.

- omega_values:

  Numeric matrix. Optional matrix of starting values for `omega`.

- omega_free:

  Logical matrix. Optional matrix of free (`TRUE`) parameters for
  `omega`.

- omega_lbound:

  Numeric matrix. Optional matrix of lower bound values for `omega`.

- omega_ubound:

  Numeric matrix. Optional matrix of upper bound values for `omega`.

- psi_diag:

  Logical. If `psi_diag = TRUE`, `psi` is a diagonal matrix. If
  `psi_diag = FALSE`, `psi` is a symmetric matrix.

- psi_d_free:

  Logical vector indicating free/fixed status of the elements of
  `psi_d`. If `NULL`, all element of `psi_d` are free.

- psi_d_values:

  Numeric vector with starting values for `psi_d`. If `NULL`, defaults
  to a vector of ones.

- psi_d_lbound:

  Numeric vector with lower bounds for `psi_d`. If `NULL`, no lower
  bounds are set.

- psi_d_ubound:

  Numeric vector with upper bounds for `psi_d`. If `NULL`, no upper
  bounds are set.

- psi_l_free:

  Logical matrix indicating which strictly-lower-triangular elements of
  `psi_l` are free. Ignored if `psi_diag = TRUE`.

- psi_l_values:

  Numeric matrix of starting values for the strictly-lower-triangular
  elements of `psi_l`. If `NULL`, defaults to a null matrix.

- psi_l_lbound:

  Numeric matrix with lower bounds for `psi_l`. If `NULL`, no lower
  bounds are set.

- psi_l_ubound:

  Numeric matrix with upper bounds for `psi_l`. If `NULL`, no upper
  bounds are set.

- check_estimates:

  Logical. Check elements of `v` for positive definiteness. If the test
  fails, the function generates a near positive definite matrix to
  replace the original using
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html).

- effects:

  Logical. If `effects = TRUE`, include estimates of the dynamic effects
  matrix, if available. If `effects = FALSE`, exclude estimates of the
  dynamic effects matrix.

- set_point:

  Logical. If `set_point = TRUE`, include estimates of the set-point
  vector, if available. If `set_point = FALSE`, exclude estimates of the
  set-point vector.

- int_meas:

  Logical. If `int_meas = TRUE`, include estimates of the measurement
  intercept vector, if available. If `int_meas = FALSE`, exclude
  estimates of the measurement intercept vector.

- int_dyn:

  Logical. If `int_dyn = TRUE`, include estimates of the dynamic process
  intercept vector, if available. If `int_dyn = FALSE`, exclude
  estimates of the dynamic process intercept vector.

- cov_meas:

  Logical. If `cov_meas = TRUE`, include estimates of the measurement
  error covariance matrix, if available. If `cov_meas = FALSE`, exclude
  estimates of the measurement error covariance matrix.

- cov_dyn:

  Logical. If `cov_dyn = TRUE`, include estimates of the process noise
  covariance matrix, if available. If `cov_dyn = FALSE`, exclude
  estimates of the process noise covariance matrix.

- diag_cov:

  Character string. If `diag_cov = "var"`, `cov_dyn` and `cov_meas` are
  in the original metric variance/covariance metric. If
  `diag_cov = "logvar"`, the diagonal elements of `cov_dyn` and
  `cov_meas` are the log of the variances and the off-diagonal elements
  are the elements in `L` in the `LDL'` decomposition. If
  `diag_cov = "softplusvar"`, the diagonal elements of `cov_dyn` and
  `cov_meas` are the softplus of the variances and the off-diagonal
  elements are the elements in `L` in the `LDL'` decomposition.

- converged:

  Logical. Only include converged cases.

- vanishing_theta:

  Logical. Test for measurement error variance going to zero if
  `converged = TRUE`.

- theta_tol:

  Numeric. Tolerance for vanishing theta test if `converged` and
  `theta_tol` are `TRUE`.

- robust_v:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix in stage 1. If `FALSE`, use normal theory sampling
  variance-covariance matrix in stage 1.

- robust:

  Logical. If `TRUE`, calculate robust (sandwich) sampling
  variance-covariance matrix in stage 2.

- lb:

  Logical. If `TRUE`, calculate likelihood based confidence intervals in
  stage 2.

- alpha:

  NUmeric. Alpha for test of significance and confidence intervals.

- tries_explore:

  Integer. Number of extra tries for the wide exploration phase using
  [`OpenMx::mxTryHardWideSearch()`](https://rdrr.io/pkg/OpenMx/man/mxTryHard.html)
  with `checkHess = FALSE`.

- tries_local:

  Integer. Number of extra tries for local polishing via
  [`OpenMx::mxTryHard()`](https://rdrr.io/pkg/OpenMx/man/mxTryHard.html)
  when gradients remain above tolerance.

- max_attempts:

  Integer. Maximum number of remediation attempts after the first
  Hessian computation fails the criteria. Each attempt may nudge off
  bounds, refit locally without the Hessian, and, on the last attempt,
  relax bounds.

- grad_tol:

  Numeric. Tolerance for the maximum absolute gradient. Smaller values
  are stricter.

- hess_tol:

  Numeric. Minimum allowable Hessian eigenvalue. Smaller values are less
  strict.

- eps:

  Numeric. Proximity threshold to detect parameters on their bounds and
  to nudge them inward by `10 * eps`.

- factor:

  Numeric. Multiplicative factor to relax parameter bounds on the final
  remediation attempt. Lower bounds are divided by `factor` and upper
  bounds are multiplied by `factor`.

- abs_bnd_tol:

  Numeric scalar. Absolute tolerance used when comparing parameter
  values to bounds.

- rel_bnd_tol:

  Numeric scalar. Relative tolerance multiplier.

- silent:

  Logical. If `TRUE`, suppresses messages during the model fitting
  stage.

- seed:

  Random seed for reproducibility.

- ncores:

  Positive integer. Number of cores to use.

## References

Cheung, M. W.-L. (2015). *Meta-analysis: A structural equation modeling
approach*. Wiley.
[doi:10.1002/9781118957813](https://doi.org/10.1002/9781118957813)

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2015). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika*, *81*(2), 535–549.
[doi:10.1007/s11336-014-9435-8](https://doi.org/10.1007/s11336-014-9435-8)

## See also

Other Meta-Analysis of VAR Functions:
[`Meta()`](https://github.com/jeksterslab/metaDyn/reference/Meta.md)

## Author

Ivan Jacob Agaloos Pesigan
