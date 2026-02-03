.MetaMxMatricesParameters <- function(p,
                                      m,
                                      r,
                                      ynames,
                                      xnames,
                                      znames,
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
                                      gamma_values,
                                      gamma_free,
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
                                      psi_l_ubound) {
  alpha <- .MetaAlpha(
    p = p,
    ynames = ynames,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound
  )
  if (random) {
    tau_sqr <- .MetaTauSqr(
      p = p,
      tau_sqr_diag = tau_sqr_diag,
      tau_sqr_d_free = tau_sqr_d_free,
      tau_sqr_d_values = tau_sqr_d_values,
      tau_sqr_d_lbound = tau_sqr_d_lbound,
      tau_sqr_d_ubound = tau_sqr_d_ubound,
      tau_sqr_l_free = tau_sqr_l_free,
      tau_sqr_l_values = tau_sqr_l_values,
      tau_sqr_l_lbound = tau_sqr_l_lbound,
      tau_sqr_l_ubound = tau_sqr_l_ubound
    )
    i_sqr <- .MetaISqr(
      univariate = i_sqr_univariate
    )
  } else {
    tau_sqr <- NULL
    i_sqr <- NULL
  }
  if (covariate) {
    gamma <- .MetaGamma(
      p = p,
      m = m,
      ynames = ynames,
      xnames = xnames,
      gamma_values = gamma_values,
      gamma_free = gamma_free,
      gamma_lbound = gamma_lbound,
      gamma_ubound = gamma_ubound
    )
  } else {
    gamma <- NULL
  }
  if (distal) {
    kappa <- .MetaKappa(
      r = r,
      znames = znames,
      kappa_free = kappa_free,
      kappa_values = kappa_values,
      kappa_lbound = kappa_lbound,
      kappa_ubound = kappa_ubound
    )
    phi <- .MetaPhi(
      p = p,
      r = r,
      ynames = ynames,
      znames = znames,
      phi_values = phi_values,
      phi_free = phi_free,
      phi_lbound = phi_lbound,
      phi_ubound = phi_ubound
    )
    psi <- .MetaPsi(
      r = r,
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
    if (covariate) {
      omega <- .MetaOmega(
        m = m,
        r = r,
        xnames = xnames,
        znames = znames,
        omega_values = omega_values,
        omega_free = omega_free,
        omega_lbound = omega_lbound,
        omega_ubound = omega_ubound
      )
    } else {
      omega <- NULL
    }
  } else {
    kappa <- NULL
    phi <- NULL
    psi <- NULL
    omega <- NULL
  }
  out <- c(
    alpha,
    tau_sqr,
    i_sqr,
    gamma,
    kappa,
    phi,
    psi,
    omega
  )
  out[
    !sapply(
      X = out,
      FUN = is.null
    )
  ]
}
