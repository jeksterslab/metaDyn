.MetaCheckFree <- function(alpha_free,
                           alpha_values,
                           tau_sqr_diag,
                           tau_sqr_d_free,
                           tau_sqr_d_values,
                           tau_sqr_l_free,
                           tau_sqr_l_values,
                           gamma_free,
                           gamma_values,
                           kappa_free,
                           kappa_values,
                           phi_free,
                           phi_values,
                           omega_free,
                           omega_values,
                           psi_diag,
                           psi_d_free,
                           psi_d_values,
                           psi_l_free,
                           psi_l_values) {
  if (is.null(alpha_values) && isFALSE(is.null(alpha_free))) {
    if (any(!c(alpha_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`alpha_values`",
          "required for fixed and free parameters in",
          "`alpha`.",
          "\n"
        )
      )
    }
  }
  if (is.null(gamma_values) && isFALSE(is.null(gamma_free))) {
    if (any(!c(gamma_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`gamma_values`",
          "required for fixed and free parameters in",
          "`gamma`.",
          "\n"
        )
      )
    }
  }
  if (is.null(kappa_values) && isFALSE(is.null(kappa_free))) {
    if (any(!c(kappa_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`kappa_values`",
          "required for fixed and free parameters in",
          "`kappa`.",
          "\n"
        )
      )
    }
  }
  if (is.null(phi_values) && isFALSE(is.null(phi_free))) {
    if (any(!c(phi_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`phi_values`",
          "required for fixed and free parameters in",
          "`phi`.",
          "\n"
        )
      )
    }
  }
  if (is.null(omega_values) && isFALSE(is.null(omega_free))) {
    if (any(!c(omega_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`omega_values`",
          "required for fixed and free parameters in",
          "`omega`.",
          "\n"
        )
      )
    }
  }
  if (is.null(tau_sqr_d_values) && isFALSE(is.null(tau_sqr_d_free))) {
    if (any(!c(tau_sqr_d_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`tau_sqr_d_values`",
          "required for fixed and free parameters in",
          "`tau_sqr_d`.",
          "\n"
        )
      )
    }
  }
  if (isFALSE(tau_sqr_diag)) {
    if (is.null(tau_sqr_l_values) && isFALSE(is.null(tau_sqr_l_free))) {
      if (any(!c(tau_sqr_l_free))) {
        stop(
          paste(
            "\nFixed and starting values",
            "`tau_sqr_l_values`",
            "required for fixed and free parameters in",
            "`tau_sqr_l`.",
            "\n"
          )
        )
      }
    }
  }
  if (is.null(psi_d_values) && isFALSE(is.null(psi_d_free))) {
    if (any(!c(psi_d_free))) {
      stop(
        paste(
          "\nFixed and starting values",
          "`psi_d_values`",
          "required for fixed and free parameters in",
          "`psi_d`.",
          "\n"
        )
      )
    }
  }
  if (isFALSE(psi_diag)) {
    if (is.null(psi_l_values) && isFALSE(is.null(psi_l_free))) {
      if (any(!c(psi_l_free))) {
        stop(
          paste(
            "\nFixed and starting values",
            "`psi_l_values`",
            "required for fixed and free parameters in",
            "`psi_l`.",
            "\n"
          )
        )
      }
    }
  }
}
