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
  gamma_free = NULL,
  gamma_values = NULL,
  gamma_lbound = NULL,
  gamma_ubound = NULL,
  kappa_free = NULL,
  kappa_values = NULL,
  kappa_lbound = NULL,
  kappa_ubound = NULL,
  phi_free = NULL,
  phi_values = NULL,
  phi_lbound = NULL,
  phi_ubound = NULL,
  omega_free = NULL,
  omega_values = NULL,
  omega_lbound = NULL,
  omega_ubound = NULL,
  psi_diag = FALSE,
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
  set_point = TRUE,
  int_meas = TRUE,
  int_dyn = TRUE,
  cov_meas = TRUE,
  cov_dyn = TRUE,
  robust_v = FALSE,
  robust = FALSE,
  alpha = 0.05,
  seed = NULL,
  tries_explore = 100,
  tries_local = 100,
  max_attempts = 10,
  silent = FALSE,
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

- phi_free:

  Logical matrix. Optional matrix of free (`TRUE`) parameters for `phi`.

- phi_values:

  Numeric matrix. Optional matrix of starting values for `phi`.

- phi_lbound:

  Numeric matrix. Optional matrix of lower bound values for `phi`.

- phi_ubound:

  Numeric matrix. Optional matrix of upper bound values for `phi`.

- omega_free:

  Logical matrix. Optional matrix of free (`TRUE`) parameters for
  `omega`.

- omega_values:

  Numeric matrix. Optional matrix of starting values for `omega`.

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

- robust_v:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix in stage 1. If `FALSE`, use normal theory sampling
  variance-covariance matrix in stage 1.

- robust:

  Logical. If `TRUE`, calculate robust (sandwich) sampling
  variance-covariance matrix in stage 2.

- alpha:

  NUmeric. Alpha for test of significance and confidence intervals.

- seed:

  Random seed for reproducibility.

- tries_explore:

  Integer. Number of extra tries for the wide exploration phase.

- tries_local:

  Integer. Number of extra tries for local polishing.

- max_attempts:

  Integer. Maximum number of remediation attempts after the first
  Hessian computation fails the criteria.

- silent:

  Logical. If `TRUE`, suppresses messages during the model fitting
  stage.

- ncores:

  Positive integer. Number of cores to use.

## Value

Returns an object of class `metadynmeta` which is a list with the
following elements:

- call:

  Function call.

- args:

  List of function arguments.

- fun:

  Function used ("Meta").

- output:

  A fitted OpenMx model.

- robust:

  Output from
  [`OpenMx::imxRobustSE()`](https://rdrr.io/pkg/OpenMx/man/imxRobustSE.html)
  with argument `details = TRUE` if `robust = TRUE`.

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

## Examples

``` r
# \donttest{
if (requireNamespace("simStateSpace")) {
  # Generate data using the simStateSpace package-------------------------
  library(simStateSpace)
  set.seed(42)
  n <- 5
  time <- 100
  p <- 2
  alpha <- rep(x = 0, times = p)
  beta <- 0.50 * diag(p)
  psi <- 0.001 * diag(p)
  psi_l <- t(chol(psi))
  mu0 <- SSMMeanEta(
    beta = beta,
    alpha = alpha
  )
  sigma0 <- SSMCovEta(
    beta = beta,
    psi = psi
  )
  sigma0_l <- t(chol(sigma0))
  sim <- SimSSMVARFixed(
    n = n,
    time = time,
    mu0 = mu0,
    sigma0_l = sigma0_l,
    alpha = alpha,
    beta = beta,
    psi_l = psi_l
  )
  data <- as.data.frame(sim)

  # Stage 1---------------------------------------------------------------
  library(fitVARMxID)
  stage1 <- FitVARMxID(
    data = data,
    observed = paste0("y", seq_len(p)),
    id = "id",
    center = TRUE
  )
  summary(stage1)
  # Stage 2---------------------------------------------------------------
  # Meta-analyze set point vector and matrix of lagged-effects
  library(metaDyn)
  stage2 <- MetaVARMx(
    object = stage1,
    random = FALSE,
    effects = TRUE,
    set_point = TRUE,
    int_meas = FALSE,
    int_dyn = FALSE,
    cov_meas = FALSE,
    cov_dyn = FALSE
  )
  # Methods for the output of the MetaVARMx() function
  print(stage2)
  summary(stage2)
  coef(stage2)
  vcov(stage2)
  confint(stage2)
  extract(stage2, what = "alpha")
}
#> Running DTVAR_ID1 with 9 parameters
#> 
#> Beginning initial fit attempt
#> Running DTVAR_ID1 with 9 parameters
#> 
#>  Lowest minimum so far:  -822.060942313193
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-822.06094 (started at 367.82482)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> 0.331612745243841,-0.0197631400065759,0.0270181292624815,0.539260086092713,0.0057946951525851,0.007106199232002,0.0623713422440293,-6.91195717901486,-6.98806630902827
#> Running DTVAR_ID2 with 9 parameters
#> 
#> Beginning initial fit attempt
#> Running DTVAR_ID2 with 9 parameters
#> 
#>  Lowest minimum so far:  -834.035953636412
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-834.03595 (started at 367.84346)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> 0.546679030340841,-0.0201020488638611,-0.0487574618599386,0.297450290600703,-0.00284153166188121,-0.0199122101192943,0.176134255014816,-7.14253861726531,-6.87719884724769
#> Running DTVAR_ID3 with 9 parameters
#> 
#> Beginning initial fit attempt
#> Running DTVAR_ID3 with 9 parameters
#> 
#>  Lowest minimum so far:  -808.550007472268
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-808.55001 (started at 367.87017)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> 0.528105561142997,-0.332501283495961,0.0108536161301484,0.361419705300276,-0.00426229999671594,-0.00173173484843314,-0.21387789357269,-6.9657202697529,-6.80024208411447
#> Running DTVAR_ID4 with 9 parameters
#> 
#> Beginning initial fit attempt
#> Running DTVAR_ID4 with 9 parameters
#> 
#>  Lowest minimum so far:  -845.339407101868
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-845.33941 (started at 367.80199)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> 0.160844162100514,-0.121824766994102,-0.0298919550660853,0.481815653357755,-0.0130955259827441,0.00326100822229185,0.00601848006477796,-7.23151794331632,-6.89984913049538
#> Running DTVAR_ID5 with 9 parameters
#> 
#> Beginning initial fit attempt
#> Running DTVAR_ID5 with 9 parameters
#> 
#>  Lowest minimum so far:  -832.354548096018
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-832.35455 (started at 367.7976)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> 0.299055985152669,0.0348621727088368,-0.0801687055865619,0.465890390856868,-0.00111115045695215,0.00331563195705935,-0.166419668597626,-6.91180267781828,-7.08999323366805
#> Running Model with 6 parameters
#> 
#> Beginning initial fit attempt
#> Running Model with 6 parameters
#> 
#>  Lowest minimum so far:  -93.184957280109
#> 
#> Solution found
#> 
#> 
#>  Solution found!  Final fit=-93.184957 (started at -91.256471)  (1 attempt(s): 1 valid, 0 errors)
#>  Start values from best fit:
#> -0.00428470046284236,-0.00406757081116166,0.386594914967842,-0.0851794164286438,-0.0234341474339058,0.444138933642437
#> Call:
#> MetaVARMx(object = stage1, random = FALSE, effects = TRUE, set_point = TRUE, 
#>     int_meas = FALSE, int_dyn = FALSE, cov_meas = FALSE, cov_dyn = FALSE)
#> 
#> Status code:
#> 0
#> 
#> CI type:
#> "normal"
#> 
#>                est     se       z      p    2.5%   97.5%
#> alpha[1,1] -0.0043 0.0020 -2.1572 0.0310 -0.0082 -0.0004
#> alpha[2,1] -0.0041 0.0025 -1.6581 0.0973 -0.0089  0.0007
#> alpha[3,1]  0.3866 0.0404  9.5726 0.0000  0.3074  0.4657
#> alpha[4,1] -0.0852 0.0426 -2.0016 0.0453 -0.1686 -0.0018
#> alpha[5,1] -0.0234 0.0366 -0.6403 0.5220 -0.0952  0.0483
#> alpha[6,1]  0.4441 0.0389 11.4250 0.0000  0.3679  0.5203
#>           alpha
#> y1 -0.004284700
#> y2 -0.004067571
#> y3  0.386594915
#> y4 -0.085179416
#> y5 -0.023434147
#> y6  0.444138934
# }
```
