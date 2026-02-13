#' Fit Multivariate Meta-Analysis (metaVAR via OpenMx)
#'
#' This function estimates fixed-, random-, or mixed-effects
#' meta-analytic parameters using per-individual coefficient estimates
#' and their sampling variance-covariance matrices.
#' Optionally, it fits distal-outcome models
#' in which between-person outcomes are regressed on
#' between-person covariates and the meta-analyzed parameters/effect sizes.
#' This function uses the estimated coefficients and
#' sampling variance-covariance matrix
#' from each individual fitted using the
#' [fitVARMxID::FitVARMxID()] function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param object Output of the [fitVARMxID::FitVARMxID()] function.
#' @param effects Logical.
#'   If `effects = TRUE`,
#'   include estimates of the dynamic effects matrix, if available.
#'   If `effects = FALSE`,
#'   exclude estimates of the dynamic effects matrix.
#' @param set_point Logical.
#'   If `set_point = TRUE`,
#'   include estimates of the set-point vector, if available.
#'   If `set_point = FALSE`,
#'   exclude estimates of the set-point vector.
#' @param int_meas Logical.
#'   If `int_meas = TRUE`,
#'   include estimates of the measurement intercept vector, if available.
#'   If `int_meas = FALSE`,
#'   exclude estimates of the measurement intercept vector.
#' @param int_dyn Logical.
#'   If `int_dyn = TRUE`,
#'   include estimates of the dynamic process intercept vector, if available.
#'   If `int_dyn = FALSE`,
#'   exclude estimates of the dynamic process intercept vector.
#' @param cov_meas Logical.
#'   If `cov_meas = TRUE`,
#'   include estimates of the measurement error covariance matrix, if available.
#'   If `cov_meas = FALSE`,
#'   exclude estimates of the measurement error covariance matrix.
#' @param cov_dyn Logical.
#'   If `cov_dyn = TRUE`,
#'   include estimates of the process noise covariance matrix, if available.
#'   If `cov_dyn = FALSE`,
#'   exclude estimates of the process noise covariance matrix.
#' @param diag_cov Character string.
#'   If `diag_cov = "var"`,
#'   `cov_dyn` and `cov_meas`
#'   are in the original metric variance/covariance metric.
#'   If `diag_cov = "logvar"`,
#'   the diagonal elements of `cov_dyn` and `cov_meas`
#'   are the log of the variances
#'   and the off-diagonal elements are the elements in `L`
#'   in the `LDL'` decomposition.
#'   If `diag_cov = "softplusvar"`,
#'   the diagonal elements of `cov_dyn` and `cov_meas`
#'   are the softplus of the variances
#'   and the off-diagonal elements are the elements in `L`
#'   in the `LDL'` decomposition.
#' @param converged Logical.
#'   Only include converged cases.
#' @param vanishing_theta Logical.
#'   Test for measurement error variance going to zero
#'   if `converged = TRUE`.
#' @param theta_tol Numeric.
#'   Tolerance for vanishing theta test
#'   if `converged` and `theta_tol` are `TRUE`.
#' @param robust_v Logical.
#'   If `TRUE`,
#'   use robust (sandwich) sampling variance-covariance matrix
#'   in stage 1.
#'   If `FALSE`,
#'   use normal theory sampling variance-covariance matrix
#'   in stage 1.
#' @inheritParams Meta
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
#' @export
MetaVARMx <- function(object,
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
                      grad_tol = 1e-2,
                      hess_tol = 1e-8,
                      eps = 1e-6,
                      factor = 10,
                      abs_bnd_tol = 1e-6,
                      rel_bnd_tol = 1e-4,
                      silent = FALSE,
                      seed = NULL,
                      ncores = NULL) {
  stopifnot(
    inherits(
      object,
      "varmxid"
    )
  )
  if (converged) {
    fit_converged <- fitVARMxID:::converged.varmxid(
      object = object,
      grad_tol = grad_tol,
      hess_tol = hess_tol,
      vanishing_theta = vanishing_theta,
      theta_tol = theta_tol,
      prop = FALSE
    )
    object$output <- object$output[
      which(
        fit_converged
      )
    ]
  }
  y <- fitVARMxID:::coef.varmxid(
    object = object,
    mu = set_point,
    alpha = int_dyn,
    beta = effects,
    nu = int_meas,
    psi = cov_dyn,
    theta = cov_meas,
    diag_cov = diag_cov,
    converged = FALSE
  )
  v <- fitVARMxID:::vcov.varmxid(
    object = object,
    mu = set_point,
    alpha = int_dyn,
    beta = effects,
    nu = int_meas,
    psi = cov_dyn,
    theta = cov_meas,
    diag_cov = diag_cov,
    converged = FALSE,
    robust = robust_v
  )
  if (!is.null(x) && converged) {
    x <- x[
      fit_converged
    ]
  }
  if (!is.null(z) && converged) {
    z <- z[
      fit_converged
    ]
  }
  out <- Meta(
    y = y,
    v = v,
    x = x,
    z = z,
    random = random,
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
    lb = lb,
    alpha = alpha,
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
  out$call <- match.call()
  out$fun <- "MetaVARMx"
  out
}
