#' Fit Multivariate Meta-Analysis
#'
#' This function estimates fixed-, random-, or mixed-effects
#' meta-analytic parameters using per-individual coefficient estimates
#' and their sampling variance-covariance matrices.
#' Optionally, it fits distal-outcome models
#' in which between-person outcomes are regressed on
#' between-person covariates and the meta-analyzed parameters/effect sizes.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param y A list.
#'   Each element of the list is a numeric vector
#'   of estimated coefficients.
#' @param v A list.
#'   Each element of the list
#'   is a sampling variance-covariance matrix of `y`.
#' @param x An optional list.
#'   Each element of the list is a numeric vector of covariates.
#' @param z An optional list.
#'   Each element of the list is a numeric vector of distal outcomes.
#' @param random Logical.
#'   If `random = TRUE`,
#'   estimates random effects.
#'   If `random = FALSE`,
#'   `tau_sqr` is a null matrix.
#' @param alpha_free Logical vector.
#'   Optional vector of free (`TRUE`) parameters for `alpha`.
#' @param alpha_values Numeric vector.
#'   Optional vector of starting values for `alpha`.
#' @param alpha_lbound Numeric vector.
#'   Optional vector of lower bound values for `alpha`.
#' @param alpha_ubound Numeric vector.
#'   Optional vector of upper bound values for `alpha`.
#' @param tau_sqr_diag Logical.
#'   If `tau_sqr_diag = TRUE`,
#'   `tau_sqr` is a diagonal matrix.
#'   If `tau_sqr_diag = FALSE`,
#'   `tau_sqr` is a symmetric matrix.
#' @param tau_sqr_d_free Logical vector
#'   indicating free/fixed status of the elements of `tau_sqr_d`.
#'   If `NULL`, all element of `tau_sqr_d` are free.
#' @param tau_sqr_d_values Numeric vector
#'   with starting values for `tau_sqr_d`.
#'   If `NULL`, defaults to a vector of ones.
#' @param tau_sqr_d_lbound Numeric vector
#'   with lower bounds for `tau_sqr_d`.
#'   If `NULL`, no lower bounds are set.
#' @param tau_sqr_d_ubound Numeric vector
#'   with upper bounds for `tau_sqr_d`.
#'   If `NULL`, no upper bounds are set.
#' @param tau_sqr_l_free Logical matrix
#'   indicating which strictly-lower-triangular elements
#'   of `tau_sqr_l` are free.
#'   Ignored if `tau_sqr_diag = TRUE`.
#' @param tau_sqr_l_values Numeric matrix
#'   of starting values for the strictly-lower-triangular elements
#'   of `tau_sqr_l`.
#'   If `NULL`, defaults to a null matrix.
#' @param tau_sqr_l_lbound Numeric matrix
#'   with lower bounds for `tau_sqr_l`.
#'   If `NULL`, no lower bounds are set.
#' @param tau_sqr_l_ubound Numeric matrix
#'   with upper bounds for `tau_sqr_l`.
#'   If `NULL`, no upper bounds are set.
#' @param i_sqr_univariate Logical.
#'   If `i_sqr_univariate = TRUE`, use the univariate formula for \eqn{I^2}.
#'   If `i_sqr_univariate = FALSE`, use the multivariate formula for \eqn{I^2}.
#' @param gamma_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `gamma`.
#' @param gamma_values Numeric matrix.
#'   Optional matrix of starting values for `gamma`.
#' @param gamma_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `gamma`.
#' @param gamma_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `gamma`.
#' @param kappa_free Logical vector.
#'   Optional vector of free (`TRUE`) parameters for `kappa`.
#' @param kappa_values Numeric vector.
#'   Optional vector of starting values for `kappa`.
#' @param kappa_lbound Numeric vector.
#'   Optional vector of lower bound values for `kappa`.
#' @param kappa_ubound Numeric vector.
#'   Optional vector of upper bound values for `kappa`.
#' @param phi_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `phi`.
#' @param phi_values Numeric matrix.
#'   Optional matrix of starting values for `phi`.
#' @param phi_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `phi`.
#' @param phi_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `phi`.
#' @param omega_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `omega`.
#' @param omega_values Numeric matrix.
#'   Optional matrix of starting values for `omega`.
#' @param omega_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `omega`.
#' @param omega_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `omega`.
#' @param psi_diag Logical.
#'   If `psi_diag = TRUE`,
#'   `psi` is a diagonal matrix.
#'   If `psi_diag = FALSE`,
#'   `psi` is a symmetric matrix.
#' @param psi_d_free Logical vector
#'   indicating free/fixed status of the elements of `psi_d`.
#'   If `NULL`, all element of `psi_d` are free.
#' @param psi_d_values Numeric vector
#'   with starting values for `psi_d`.
#'   If `NULL`, defaults to a vector of ones.
#' @param psi_d_lbound Numeric vector
#'   with lower bounds for `psi_d`.
#'   If `NULL`, no lower bounds are set.
#' @param psi_d_ubound Numeric vector
#'   with upper bounds for `psi_d`.
#'   If `NULL`, no upper bounds are set.
#' @param psi_l_free Logical matrix
#'   indicating which strictly-lower-triangular elements
#'   of `psi_l` are free.
#'   Ignored if `psi_diag = TRUE`.
#' @param psi_l_values Numeric matrix
#'   of starting values for the strictly-lower-triangular elements
#'   of `psi_l`.
#'   If `NULL`, defaults to a null matrix.
#' @param psi_l_lbound Numeric matrix
#'   with lower bounds for `psi_l`.
#'   If `NULL`, no lower bounds are set.
#' @param psi_l_ubound Numeric matrix
#'   with upper bounds for `psi_l`.
#'   If `NULL`, no upper bounds are set.
#' @param check_estimates Logical.
#'   Check elements of `v` for positive definiteness.
#'   If the test fails, the function generates a near positive definite matrix
#'   to replace the original using [Matrix::nearPD()].
#' @param robust Logical.
#'   If `TRUE`,
#'   use robust (sandwich) sampling variance-covariance matrix
#'   in stage 2.
#'   If `FALSE`,
#'   use normal theory sampling variance-covariance matrix
#'   in stage 2.
#' @param tries_explore Integer.
#'   Number of extra tries for the wide exploration
#'   phase using `OpenMx::mxTryHardWideSearch()` with `checkHess = FALSE`.
#' @param tries_local Integer.
#'   Number of extra tries for local polishing via
#'   `OpenMx::mxTryHard()` when gradients remain above tolerance.
#' @param max_attempts Integer.
#'   Maximum number of remediation attempts
#'   after the first Hessian computation fails the criteria.
#'   Each attempt may nudge off bounds,
#'   refit locally without the Hessian, and,
#'   on the last attempt, relax bounds.
#' @param grad_tol Numeric.
#'   Tolerance for the maximum absolute gradient.
#'   Smaller values are stricter.
#' @param hess_tol Numeric.
#'   Minimum allowable Hessian eigenvalue.
#'   Smaller values are less strict.
#' @param eps Numeric.
#'   Proximity threshold to detect parameters on their bounds
#'   and to nudge them inward by `10 * eps`.
#' @param factor Numeric.
#'   Multiplicative factor to relax parameter bounds
#'   on the final remediation attempt.
#'   Lower bounds are divided by `factor` and
#'   upper bounds are multiplied by `factor`.
#' @param abs_bnd_tol Numeric scalar. Absolute tolerance used when comparing
#'   parameter values to bounds.
#' @param rel_bnd_tol Numeric scalar. Relative tolerance multiplier.
#' @param seed Random seed for reproducibility.
#' @param silent Logical.
#'   If `TRUE`, suppresses messages during the model fitting stage.
#' @param ncores Positive integer.
#'   Number of cores to use.
#'
#' @return Returns an object of class `metadynmeta` which is
#'   a list with the following elements:
#'   \describe{
#'     \item{call}{Function call.}
#'     \item{args}{List of function arguments.}
#'     \item{fun}{Function used ("Meta").}
#'     \item{output}{A fitted OpenMx model.}
#'     \item{robust}{Output from [OpenMx::imxRobustSE()]
#'         with argument `details = TRUE`.}
#'   }
#'
#' @references
#' Cheung, M. W.-L. (2015).
#' *Meta-analysis: A structural equation modeling approach*.
#' Wiley.
#' \doi{10.1002/9781118957813}
#'
#' Neale, M. C., Hunter, M. D., Pritikin, J. N.,
#' Zahery, M., Brick, T. R., Kirkpatrick, R. M., Estabrook, R.,
#' Bates, T. C., Maes, H. H., & Boker, S. M. (2015).
#' OpenMx 2.0: Extended structural equation and statistical modeling.
#' *Psychometrika*,
#' *81*(2), 535–549.
#' \doi{10.1007/s11336-014-9435-8}
#'
#' @family Meta-Analysis of VAR Functions
#' @keywords metaDyn meta
#' @import OpenMx
#' @importFrom stats coef vcov
#' @export
Meta <- function(y,
                 v,
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
                 robust = FALSE,
                 tries_explore = 100,
                 tries_local = 100,
                 max_attempts = 10,
                 grad_tol = 1e-2,
                 hess_tol = 1e-8,
                 eps = 1e-6,
                 factor = 10,
                 abs_bnd_tol = 1e-6,
                 rel_bnd_tol = 1e-4,
                 silent = FALSE,
                 seed = NULL,
                 ncores = NULL) {
  p <- length(y[[1]])
  n <- length(y)
  stopifnot(
    length(v) == length(y)
  )
  if (is.null(x)) {
    m <- NULL
    covariate <- FALSE
  } else {
    m <- length(x[[1]])
    stopifnot(
      length(x) == length(y)
    )
    covariate <- TRUE
  }
  if (is.null(z)) {
    r <- NULL
    distal <- FALSE
  } else {
    r <- length(z[[1]])
    stopifnot(
      length(z) == length(y)
    )
    distal <- TRUE
  }
  if (check_estimates) {
    estimates <- .CheckEstimates(
      y = y,
      v = v,
      p = p,
      ynames = paste0("y", seq_len(p)),
      ncores = ncores
    )
    y <- estimates$y
    v <- estimates$v
  }
  args <- list(
    y = y,
    v = v,
    x = x,
    z = z,
    p = p,
    m = m,
    r = r,
    n = n,
    random = random,
    covariate = covariate,
    distal = distal,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound,
    tau_sqr_diag = tau_sqr_diag,
    tau_sqr_d_free = tau_sqr_d_free,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_d_lbound = tau_sqr_d_lbound,
    tau_sqr_d_ubound = tau_sqr_d_ubound,
    tau_sqr_l_free = tau_sqr_l_free,
    tau_sqr_l_values = tau_sqr_l_values,
    tau_sqr_l_lbound = tau_sqr_l_lbound,
    tau_sqr_l_ubound = tau_sqr_l_ubound,
    i_sqr_univariate = i_sqr_univariate,
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    gamma_lbound = gamma_lbound,
    gamma_ubound = gamma_ubound,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    kappa_lbound = kappa_lbound,
    kappa_ubound = kappa_ubound,
    phi_values = phi_values,
    phi_free = phi_free,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_values = omega_values,
    omega_free = omega_free,
    omega_lbound = omega_lbound,
    omega_ubound = omega_ubound,
    psi_diag = psi_diag,
    psi_d_free = psi_d_free,
    psi_d_values = psi_d_values,
    psi_d_lbound = psi_d_lbound,
    psi_d_ubound = psi_d_ubound,
    psi_l_free = psi_l_free,
    psi_l_values = psi_l_values,
    psi_l_lbound = psi_l_lbound,
    psi_l_ubound = psi_l_ubound,
    check_estimates = check_estimates,
    robust = robust,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    grad_tol = grad_tol,
    hess_tol = hess_tol,
    eps = eps,
    factor = factor,
    abs_bnd_tol = abs_bnd_tol,
    rel_bnd_tol = rel_bnd_tol,
    seed = seed,
    silent = silent,
    ncores = ncores
  )
  output <- .MetaFit(
    y = y,
    v = v,
    x = x,
    z = z,
    p = p,
    m = m,
    r = r,
    n = n,
    random = random,
    covariate = covariate,
    distal = distal,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound,
    tau_sqr_diag = tau_sqr_diag,
    tau_sqr_d_free = tau_sqr_d_free,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_d_lbound = tau_sqr_d_lbound,
    tau_sqr_d_ubound = tau_sqr_d_ubound,
    tau_sqr_l_free = tau_sqr_l_free,
    tau_sqr_l_values = tau_sqr_l_values,
    tau_sqr_l_lbound = tau_sqr_l_lbound,
    tau_sqr_l_ubound = tau_sqr_l_ubound,
    i_sqr_univariate = i_sqr_univariate,
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    gamma_lbound = gamma_lbound,
    gamma_ubound = gamma_ubound,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    kappa_lbound = kappa_lbound,
    kappa_ubound = kappa_ubound,
    phi_values = phi_values,
    phi_free = phi_free,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_values = omega_values,
    omega_free = omega_free,
    omega_lbound = omega_lbound,
    omega_ubound = omega_ubound,
    psi_diag = psi_diag,
    psi_d_free = psi_d_free,
    psi_d_values = psi_d_values,
    psi_d_lbound = psi_d_lbound,
    psi_d_ubound = psi_d_ubound,
    psi_l_free = psi_l_free,
    psi_l_values = psi_l_values,
    psi_l_lbound = psi_l_lbound,
    psi_l_ubound = psi_l_ubound,
    intervals = TRUE,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    grad_tol = grad_tol,
    hess_tol = hess_tol,
    eps = eps,
    factor = factor,
    abs_bnd_tol = abs_bnd_tol,
    rel_bnd_tol = rel_bnd_tol,
    seed = seed,
    silent = silent,
    ncores = ncores
  )
  if (robust) {
    utils::capture.output(
      suppressMessages(
        suppressWarnings(
          sandwich <- OpenMx::imxRobustSE(
            model = output,
            details = TRUE
          )
        )
      )
    )
  } else {
    sandwich <- NULL
  }
  out <- list(
    call = match.call(),
    args = args,
    fun = "Meta",
    output = output,
    robust = sandwich
  )
  class(out) <- c(
    "metadynmeta",
    class(out)
  )
  out
}
