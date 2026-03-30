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
                                      gamma_free,
                                      gamma_values,
                                      gamma_lbound,
                                      gamma_ubound,
                                      kappa_free,
                                      kappa_values,
                                      kappa_lbound,
                                      kappa_ubound,
                                      phi_free,
                                      phi_values,
                                      phi_lbound,
                                      phi_ubound,
                                      omega_free,
                                      omega_values,
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
  mat_alpha <- .MetaAlpha(
    p = p,
    ynames = ynames,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound
  )
  if (random) {
    mat_tau_sqr <- .MetaTauSqr(
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
    mat_i_sqr <- .MetaISqr(
      p = p
    )
  } else {
    mat_tau_sqr <- NULL
    mat_i_sqr <- NULL
  }
  if (covariate) {
    mat_gamma <- .MetaGamma(
      p = p,
      m = m,
      ynames = ynames,
      xnames = xnames,
      gamma_free = gamma_free,
      gamma_values = gamma_values,
      gamma_lbound = gamma_lbound,
      gamma_ubound = gamma_ubound
    )
  } else {
    mat_gamma <- NULL
  }
  if (distal) {
    mat_kappa <- .MetaKappa(
      r = r,
      znames = znames,
      kappa_free = kappa_free,
      kappa_values = kappa_values,
      kappa_lbound = kappa_lbound,
      kappa_ubound = kappa_ubound
    )
    mat_phi <- .MetaPhi(
      p = p,
      r = r,
      ynames = ynames,
      znames = znames,
      phi_free = phi_free,
      phi_values = phi_values,
      phi_lbound = phi_lbound,
      phi_ubound = phi_ubound
    )
    mat_psi <- .MetaPsi(
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
      mat_omega <- .MetaOmega(
        m = m,
        r = r,
        xnames = xnames,
        znames = znames,
        omega_free = omega_free,
        omega_values = omega_values,
        omega_lbound = omega_lbound,
        omega_ubound = omega_ubound
      )
    } else {
      mat_omega <- NULL
    }
  } else {
    mat_kappa <- NULL
    mat_phi <- NULL
    mat_psi <- NULL
    mat_omega <- NULL
  }
  if (covariate && distal) {
    mat_indirect <- .MetaIndirect(
      ynames = ynames,
      xnames = xnames,
      znames = znames
    )
  } else {
    mat_indirect <- NULL
  }
  out <- c(
    mat_alpha,
    mat_tau_sqr,
    mat_i_sqr,
    mat_gamma,
    mat_kappa,
    mat_phi,
    mat_psi,
    mat_omega,
    mat_indirect
  )
  out[
    !sapply(
      X = out,
      FUN = is.null
    )
  ]
}
