.MetaFit <- function(y,
                     v,
                     x,
                     z,
                     p,
                     m,
                     r,
                     n,
                     random,
                     covariate,
                     distal,
                     alpha_free,
                     alpha_values,
                     alpha_lbound,
                     alpha_ubound,
                     tau_sqr_diag,
                     tau_sqr_d_free,
                     tau_sqr_d_values,
                     tau_sqr_d_lbound,
                     tau_sqr_d_ubound,
                     tau_sqr_l_free,
                     tau_sqr_l_values,
                     tau_sqr_l_lbound,
                     tau_sqr_l_ubound,
                     i_sqr_univariate,
                     gamma_free,
                     gamma_values,
                     gamma_lbound,
                     gamma_ubound,
                     kappa_free,
                     kappa_values,
                     kappa_lbound,
                     kappa_ubound,
                     phi_values,
                     phi_free,
                     phi_lbound,
                     phi_ubound,
                     omega_values,
                     omega_free,
                     omega_lbound,
                     omega_ubound,
                     psi_diag,
                     psi_d_free,
                     psi_d_values,
                     psi_d_lbound,
                     psi_d_ubound,
                     psi_l_free,
                     psi_l_values,
                     psi_l_lbound,
                     psi_l_ubound,
                     intervals,
                     tries_explore,
                     tries_local,
                     max_attempts,
                     grad_tol,
                     hess_tol,
                     eps,
                     factor,
                     abs_bnd_tol,
                     rel_bnd_tol,
                     seed,
                     silent,
                     ncores) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  threads <- OpenMx::mxOption(
    key = "Number of Threads"
  )
  on.exit(
    OpenMx::mxOption(
      key = "Number of Threads",
      value = threads
    ),
    add = TRUE
  )
  # nocov start
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (ncores > 1) {
      OpenMx::mxOption(
        key = "Number of Threads",
        value = ncores
      )
    }
  }
  # nocov end
  model <- .MetaModel(
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
    psi_l_ubound = psi_l_ubound
  )
  model <- OpenMx::mxTryHard(
    model = model,
    intervals = intervals,
    silent = silent
  )
  .MxHelperEnsureGoodHessian(
    model = model,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    grad_tol = grad_tol,
    hess_tol = hess_tol,
    eps = eps,
    factor = factor,
    abs_bnd_tol = abs_bnd_tol,
    rel_bnd_tol = rel_bnd_tol,
    silent = silent
  )
}
