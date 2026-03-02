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
#' @param robust_v Logical.
#'   If `TRUE`,
#'   use robust (sandwich) sampling variance-covariance matrix
#'   in stage 1.
#'   If `FALSE`,
#'   use normal theory sampling variance-covariance matrix
#'   in stage 1.
#' @inheritParams Meta
#' @inherit return Meta
#' @inherit references Meta
#'
#' @examples
#' \donttest{
#' # Generate data using the simStateSpace package-------------------------
#' library(simStateSpace)
#' set.seed(42)
#' n <- 5
#' time <- 100
#' p <- 2
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- 0.001 * diag(p)
#' psi_l <- t(chol(psi))
#' mu0 <- simStateSpace::SSMMeanEta(
#'   beta = beta,
#'   alpha = alpha
#' )
#' sigma0 <- simStateSpace::SSMCovEta(
#'   beta = beta,
#'   psi = psi
#' )
#' sigma0_l <- t(chol(sigma0))
#' sim <- SimSSMVARFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l
#' )
#' data <- as.data.frame(sim)
#'
#' # Stage 1---------------------------------------------------------------
#' library(fitVARMxID)
#' stage1 <- FitVARMxID(
#'   data = data,
#'   observed = paste0("y", seq_len(p)),
#'   id = "id",
#'   center = TRUE
#' )
#' summary(stage1)
#' # Stage 2---------------------------------------------------------------
#' # Meta-analyze set point vector and matrix of lagged-effects
#' stage2 <- MetaVARMx(
#'   object = stage1,
#'   random = FALSE,
#'   effects = TRUE,
#'   set_point = TRUE,
#'   int_meas = FALSE,
#'   int_dyn = FALSE,
#'   cov_meas = FALSE,
#'   cov_dyn = FALSE
#' )
#' summary(stage2)
#' }
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
                      ncores = NULL) {
  stopifnot(
    inherits(
      object,
      "varmxid"
    )
  )
  y <- fitVARMxID:::coef.varmxid(
    object = object,
    mu = set_point,
    alpha = int_dyn,
    beta = effects,
    nu = int_meas,
    psi = cov_dyn,
    theta = cov_meas,
    ncores = ncores
  )
  v <- fitVARMxID:::vcov.varmxid(
    object = object,
    mu = set_point,
    alpha = int_dyn,
    beta = effects,
    nu = int_meas,
    psi = cov_dyn,
    theta = cov_meas,
    robust = robust_v,
    ncores = ncores
  )
  out <- Meta(
    y = y,
    v = v,
    x = x[
      object$converged
    ],
    z = z[
      object$converged
    ],
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
    phi_free = phi_free,
    phi_values = phi_values,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_free = omega_free,
    omega_values = omega_values,
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
    alpha = alpha,
    seed = seed,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    silent = silent,
    ncores = ncores
  )
  out$call <- match.call()
  out$fun <- "MetaVARMx"
  out
}
