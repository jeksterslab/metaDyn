.MetaStarts <- function(raw_data,
                        fixed_x,
                        xnames,
                        ynames,
                        znames,
                        alpha_values,
                        mu_x_values,
                        sigma_x_d_values,
                        sigma_x_l_values,
                        gamma_values,
                        kappa_values,
                        phi_values,
                        omega_values,
                        tau_sqr_d_values,
                        tau_sqr_l_values,
                        psi_d_values,
                        psi_l_values) {
  .CovStartValues <- function(sigma,
                              min_diag = 0.001) {
    if (is.null(dim(sigma))) {
      sigma <- matrix(
        data = sigma,
        nrow = 1,
        ncol = 1
      )
    }
    sigma_diag <- diag(sigma)
    sigma_diag <- ifelse(
      test = is.na(sigma_diag) | sigma_diag <= min_diag,
      yes = min_diag,
      no = sigma_diag
    )
    diag(sigma) <- sigma_diag
    if (ncol(sigma) == 1) {
      out <- list(
        d_values = .MxHelperInvSoftplus(c(sigma)),
        l_values = NULL
      )
    } else {
      ldl_sigma <- .MxHelperLDL(
        x = sigma
      )
      out <- list(
        d_values = ldl_sigma$uc_d,
        l_values = ldl_sigma$s_l
      )
    }
    out
  }
  complete <- stats::complete.cases(raw_data)
  if (!any(complete)) {
    message(
      paste(
        "No complete rows.",
        "Ordinary least squares starting values cannot be provided."
      )
    )
  } else {
    raw_data_complete <- as.matrix(
      raw_data[
        complete,
        ,
        drop = FALSE
      ]
    )
    y_obs <- raw_data_complete[
      ,
      ynames,
      drop = FALSE
    ]
    if (!is.null(xnames)) {
      x_obs <- raw_data_complete[
        ,
        xnames,
        drop = FALSE
      ]
    } else {
      x_obs <- NULL
    }
    # Stochastic x starts.
    # These are independent of whether alpha/gamma starts were supplied.
    if (!fixed_x && !is.null(xnames)) {
      if (is.null(mu_x_values)) {
        mu_x_values <- colMeans(
          x_obs
        )
      }
      if (is.null(sigma_x_d_values) || is.null(sigma_x_l_values)) {
        sigma_x <- stats::var(
          x = x_obs
        )
        sigma_x_starts <- .CovStartValues(
          sigma = sigma_x
        )
        if (is.null(sigma_x_d_values)) {
          sigma_x_d_values <- sigma_x_starts$d_values
        }
        if (is.null(sigma_x_l_values)) {
          sigma_x_l_values <- sigma_x_starts$l_values
        }
      }
    }
    # y on x, or y intercept only when no x is available.
    if (is.null(xnames)) {
      if (is.null(alpha_values)) {
        alpha_values <- colMeans(
          y_obs
        )
      }
      if (is.null(tau_sqr_d_values) || is.null(tau_sqr_l_values)) {
        tau_sqr <- stats::var(
          x = y_obs
        )
        tau_sqr_starts <- .CovStartValues(
          sigma = tau_sqr
        )
        if (is.null(tau_sqr_d_values)) {
          tau_sqr_d_values <- tau_sqr_starts$d_values
        }
        if (is.null(tau_sqr_l_values)) {
          tau_sqr_l_values <- tau_sqr_starts$l_values
        }
      }
    } else {
      need_y_ols <- (
        is.null(alpha_values) ||
          is.null(gamma_values) ||
          is.null(tau_sqr_d_values) ||
          is.null(tau_sqr_l_values)
      )
      if (need_y_ols) {
        x_mat <- cbind(
          1,
          x_obs
        )
        colnames(x_mat) <- c(
          "intercept",
          xnames
        )
        qrx <- qr(
          x = x_mat
        )
        gamma_hat <- t(
          qr.coef(
            qr = qrx,
            y = y_obs
          )
        )
        if (is.null(alpha_values)) {
          alpha_values <- gamma_hat[
            ,
            "intercept",
            drop = TRUE
          ]
        }
        if (is.null(gamma_values)) {
          gamma_values <- gamma_hat[
            ,
            xnames,
            drop = FALSE
          ]
        }
        if (is.null(tau_sqr_d_values) || is.null(tau_sqr_l_values)) {
          tau_sqr <- (
            crossprod(
              y_obs - x_mat %*% t(gamma_hat)
            ) / (
              nrow(y_obs) - qrx$rank
            )
          )
          tau_sqr_starts <- .CovStartValues(
            sigma = tau_sqr
          )
          if (is.null(tau_sqr_d_values)) {
            tau_sqr_d_values <- tau_sqr_starts$d_values
          }
          if (is.null(tau_sqr_l_values)) {
            tau_sqr_l_values <- tau_sqr_starts$l_values
          }
        }
      }
    }
    # z on y, or z on x + y.
    # These are pragmatic OLS starts; the fitted model still uses
    # latent effect sizes to predict z.
    if (!is.null(znames)) {
      z_obs <- raw_data_complete[
        ,
        znames,
        drop = FALSE
      ]
      if (is.null(xnames)) {
        need_z_ols <- (
          is.null(kappa_values) ||
            is.null(phi_values) ||
            is.null(psi_d_values) ||
            is.null(psi_l_values)
        )
        if (need_z_ols) {
          y_mat <- cbind(
            1,
            y_obs
          )
          colnames(y_mat) <- c(
            "intercept",
            ynames
          )
          qry <- qr(
            x = y_mat
          )
          phi_hat <- t(
            qr.coef(
              qr = qry,
              y = z_obs
            )
          )
          if (is.null(kappa_values)) {
            kappa_values <- phi_hat[
              ,
              "intercept",
              drop = TRUE
            ]
          }
          if (is.null(phi_values)) {
            phi_values <- phi_hat[
              ,
              ynames,
              drop = FALSE
            ]
          }
          if (is.null(psi_d_values) || is.null(psi_l_values)) {
            psi <- (
              crossprod(
                z_obs - y_mat %*% t(phi_hat)
              ) / (
                nrow(z_obs) - qry$rank
              )
            )
            psi_starts <- .CovStartValues(
              sigma = psi
            )
            if (is.null(psi_d_values)) {
              psi_d_values <- psi_starts$d_values
            }
            if (is.null(psi_l_values)) {
              psi_l_values <- psi_starts$l_values
            }
          }
        }
      } else {
        need_z_ols <- (
          is.null(kappa_values) ||
            is.null(omega_values) ||
            is.null(phi_values) ||
            is.null(psi_d_values) ||
            is.null(psi_l_values)
        )
        if (need_z_ols) {
          xy_mat <- cbind(
            1,
            x_obs,
            y_obs
          )
          colnames(xy_mat) <- c(
            "intercept",
            xnames,
            ynames
          )
          qrxy <- qr(
            x = xy_mat
          )
          phi_hat <- t(
            qr.coef(
              qr = qrxy,
              y = z_obs
            )
          )
          if (is.null(kappa_values)) {
            kappa_values <- phi_hat[
              ,
              "intercept",
              drop = TRUE
            ]
          }
          if (is.null(omega_values)) {
            omega_values <- phi_hat[
              ,
              xnames,
              drop = FALSE
            ]
          }
          if (is.null(phi_values)) {
            phi_values <- phi_hat[
              ,
              ynames,
              drop = FALSE
            ]
          }
          if (is.null(psi_d_values) || is.null(psi_l_values)) {
            psi <- (
              crossprod(
                z_obs - xy_mat %*% t(phi_hat)
              ) / (
                nrow(z_obs) - qrxy$rank
              )
            )
            psi_starts <- .CovStartValues(
              sigma = psi
            )
            if (is.null(psi_d_values)) {
              psi_d_values <- psi_starts$d_values
            }
            if (is.null(psi_l_values)) {
              psi_l_values <- psi_starts$l_values
            }
          }
        }
      }
    }
  }
  list(
    alpha_values = alpha_values,
    mu_x_values = mu_x_values,
    sigma_x_d_values = sigma_x_d_values,
    sigma_x_l_values = sigma_x_l_values,
    gamma_values = gamma_values,
    kappa_values = kappa_values,
    phi_values = phi_values,
    omega_values = omega_values,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_l_values = tau_sqr_l_values,
    psi_d_values = psi_d_values,
    psi_l_values = psi_l_values
  )
}
