.MetaStarts <- function(raw_data,
                        xnames,
                        ynames,
                        znames,
                        alpha_values,
                        gamma_values,
                        kappa_values,
                        phi_values,
                        omega_values,
                        tau_sqr_d_values,
                        tau_sqr_l_values,
                        psi_d_values,
                        psi_l_values) {
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
        complete, ,
        drop = FALSE
      ]
    )
    # y on x
    if (is.null(xnames)) {
      if (is.null(alpha_values)) {
        y_mat <- raw_data_complete[
          ,
          ynames,
          drop = FALSE
        ]
        alpha_values <- colMeans(
          y_mat
        )
        tau_sqr <- stats::var(y_mat)
        if (length(ynames) == 1) {
          tau_sqr_d_values <- .MxHelperInvSoftplus(c(tau_sqr))
        } else {
          ldl_tau_sqr <- .MxHelperLDL(
            x = tau_sqr
          )
          tau_sqr_d_values <- ldl_tau_sqr$uc_d
          tau_sqr_l_values <- ldl_tau_sqr$s_l
        }
      }
    } else {
      if (is.null(alpha_values) || is.null(gamma_values)) {
        x_mat <- cbind(
          1,
          raw_data_complete[
            ,
            xnames,
            drop = FALSE
          ]
        )
        y_mat <- raw_data_complete[
          ,
          ynames,
          drop = FALSE
        ]
        qrx <- qr(
          x = x_mat
        )
        gamma_hat <- t(
          qr.coef(
            qr = qrx,
            y = y_mat
          )
        )
        alpha_values <- gamma_hat[, 1]
        gamma_values <- gamma_hat[
          ,
          xnames,
          drop = FALSE
        ]
        tau_sqr <- (
          crossprod(
            y_mat - x_mat %*% t(gamma_hat)
          ) / (
            nrow(y_mat) - qr(x_mat)$rank
          )
        )
        # check if diag is lower than epsilon
        tau_sqr_diag <- diag(tau_sqr)
        tau_sqr_diag <- ifelse(
          test = tau_sqr_diag <= 0.001,
          yes = 0.001,
          no = tau_sqr_diag
        )
        diag(tau_sqr) <- tau_sqr_diag
        if (length(ynames) == 1) {
          tau_sqr_d_values <- .MxHelperInvSoftplus(c(tau_sqr))
        } else {
          ldl_tau_sqr <- .MxHelperLDL(
            x = tau_sqr
          )
          tau_sqr_d_values <- ldl_tau_sqr$uc_d
          tau_sqr_l_values <- ldl_tau_sqr$s_l
        }
      }
    }
    # z on x + y
    if (isFALSE(is.null(znames))) {
      if (is.null(xnames)) {
        if (is.null(kappa_values) || is.null(phi_values)) {
          y_mat <- cbind(
            1,
            raw_data_complete[
              ,
              ynames,
              drop = FALSE
            ]
          )
          z_mat <- raw_data_complete[
            ,
            znames,
            drop = FALSE
          ]
          qry <- qr(
            x = y_mat
          )
          phi_hat <- t(
            qr.coef(
              qr = qry,
              y = z_mat
            )
          )
          kappa_values <- phi_hat[, 1]
          phi_values <- phi_hat[
            ,
            ynames,
            drop = FALSE
          ]
          psi <- (
            crossprod(
              z_mat - y_mat %*% t(phi_hat)
            ) / (
              nrow(z_mat) - qr(y_mat)$rank
            )
          )
          # check if diag is lower than epsilon
          psi_diag <- diag(psi)
          psi_diag <- ifelse(
            test = psi_diag <= 0.001,
            yes = 0.001,
            no = psi_diag
          )
          diag(psi) <- psi_diag
          if (length(znames) == 1) {
            psi_d_values <- .MxHelperInvSoftplus(c(psi))
          } else {
            ldl_psi <- .MxHelperLDL(
              x = psi
            )
            psi_d_values <- ldl_psi$uc_d
            psi_l_values <- ldl_psi$s_l
          }
        }
      } else {
        if (is.null(kappa_values) || is.null(phi_values)) {
          xy_mat <- cbind(
            1,
            raw_data_complete[
              ,
              c(xnames, ynames),
              drop = FALSE
            ]
          )
          z_mat <- raw_data_complete[
            ,
            znames,
            drop = FALSE
          ]
          qrxy <- qr(
            x = xy_mat
          )
          phi_hat <- t(
            qr.coef(
              qr = qrxy,
              y = raw_data_complete[
                ,
                znames,
                drop = FALSE
              ]
            )
          )
          kappa_values <- phi_hat[, 1]
          omega_values <- phi_hat[
            ,
            xnames,
            drop = FALSE
          ]
          phi_values <- phi_hat[
            ,
            ynames,
            drop = FALSE
          ]
          psi <- (
            crossprod(
              z_mat - xy_mat %*% t(phi_hat)
            ) / (
              nrow(z_mat) - qr(xy_mat)$rank
            )
          )
          ldl_psi <- .MxHelperLDL(
            x = psi
          )
          psi_d_values <- ldl_psi$uc_d
          psi_l_values <- ldl_psi$s_l
        }
      }
    }
  }
  list(
    alpha_values = alpha_values,
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
