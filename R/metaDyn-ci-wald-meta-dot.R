.CIWaldMeta <- function(object,
                        alpha) {
  mat_names <- c(
    names(object$output$algebras),
    names(object$output$matrices)
  )
  if ("alpha_vec" %in% mat_names) {
    alpha_vec <- OpenMx::mxEvalByName(
      name = "alpha_vec",
      model = object$output,
      compute = TRUE
    )
    names_alpha_vec <- rownames(alpha_vec)
    alpha_vec <- c(alpha_vec)
    names(alpha_vec) <- names_alpha_vec
    se_alpha_vec <- c(
      OpenMx::mxSE(
        x = "alpha_vec",
        model = object$output,
        silent = TRUE
      )
    )
    y0 <- .CIWald(
      est = alpha_vec,
      se = se_alpha_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    y0 <- NULL
  }
  if ("gamma_vec" %in% mat_names) {
    gamma_vec <- OpenMx::mxEvalByName(
      name = "gamma_vec",
      model = object$output,
      compute = TRUE
    )
    names_gamma_vec <- rownames(gamma_vec)
    gamma_vec <- c(gamma_vec)
    names(gamma_vec) <- names_gamma_vec
    se_gamma_vec <- c(
      OpenMx::mxSE(
        x = "gamma_vec",
        model = object$output,
        silent = TRUE
      )
    )
    y1 <- .CIWald(
      est = gamma_vec,
      se = se_gamma_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    y1 <- NULL
  }
  if ("mu_x_vec" %in% mat_names) {
    mu_x_vec <- OpenMx::mxEvalByName(
      name = "mu_x_vec",
      model = object$output,
      compute = TRUE
    )
    names_mu_x_vec <- rownames(mu_x_vec)
    mu_x_vec <- c(mu_x_vec)
    names(mu_x_vec) <- names_mu_x_vec
    se_mu_x_vec <- c(
      OpenMx::mxSE(
        x = "mu_x_vec",
        model = object$output,
        silent = TRUE
      )
    )
    x0 <- .CIWald(
      est = mu_x_vec,
      se = se_mu_x_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    x0 <- NULL
  }
  if ("beta_vec" %in% mat_names) {
    beta_vec <- OpenMx::mxEvalByName(
      name = "beta_vec",
      model = object$output,
      compute = TRUE
    )
    names_beta_vec <- rownames(beta_vec)
    beta_vec <- c(beta_vec)
    names(beta_vec) <- names_beta_vec
    se_beta_vec <- c(
      OpenMx::mxSE(
        x = "beta_vec",
        model = object$output,
        silent = TRUE
      )
    )
    yy <- .CIWald(
      est = beta_vec,
      se = se_beta_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    yy <- NULL
  }
  if ("kappa_vec" %in% mat_names) {
    kappa_vec <- OpenMx::mxEvalByName(
      name = "kappa_vec",
      model = object$output,
      compute = TRUE
    )
    names_kappa_vec <- rownames(kappa_vec)
    kappa_vec <- c(kappa_vec)
    names(kappa_vec) <- names_kappa_vec
    se_kappa_vec <- c(
      OpenMx::mxSE(
        x = "kappa_vec",
        model = object$output,
        silent = TRUE
      )
    )
    z0 <- .CIWald(
      est = kappa_vec,
      se = se_kappa_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    z0 <- NULL
  }
  if ("phi_vec" %in% mat_names) {
    phi_vec <- OpenMx::mxEvalByName(
      name = "phi_vec",
      model = object$output,
      compute = TRUE
    )
    names_phi_vec <- rownames(phi_vec)
    phi_vec <- c(phi_vec)
    names(phi_vec) <- names_phi_vec
    se_phi_vec <- c(
      OpenMx::mxSE(
        x = "phi_vec",
        model = object$output,
        silent = TRUE
      )
    )
    z1 <- .CIWald(
      est = phi_vec,
      se = se_phi_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    z1 <- NULL
  }
  if ("omega_vec" %in% mat_names) {
    omega_vec <- OpenMx::mxEvalByName(
      name = "omega_vec",
      model = object$output,
      compute = TRUE
    )
    names_omega_vec <- rownames(omega_vec)
    omega_vec <- c(omega_vec)
    names(omega_vec) <- names_omega_vec
    se_omega_vec <- c(
      OpenMx::mxSE(
        x = "omega_vec",
        model = object$output,
        silent = TRUE
      )
    )
    zx <- .CIWald(
      est = omega_vec,
      se = se_omega_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    zx <- NULL
  }
  if ("psi_vec" %in% mat_names) {
    psi_vec <- OpenMx::mxEvalByName(
      name = "psi_vec",
      model = object$output,
      compute = TRUE
    )
    names_psi_vec <- rownames(psi_vec)
    psi_vec <- c(psi_vec)
    names(psi_vec) <- names_psi_vec
    se_psi_vec <- c(
      OpenMx::mxSE(
        x = "psi_vec",
        model = object$output,
        silent = TRUE
      )
    )
    psi <- .CIWald(
      est = psi_vec,
      se = se_psi_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    psi <- NULL
  }
  if ("sigma_x_vec" %in% mat_names) {
    sigma_x_vec <- OpenMx::mxEvalByName(
      name = "sigma_x_vec",
      model = object$output,
      compute = TRUE
    )
    names_sigma_x_vec <- rownames(sigma_x_vec)
    sigma_x_vec <- c(sigma_x_vec)
    names(sigma_x_vec) <- names_sigma_x_vec
    se_sigma_x_vec <- c(
      OpenMx::mxSE(
        x = "sigma_x_vec",
        model = object$output,
        silent = TRUE
      )
    )
    sx <- .CIWald(
      est = sigma_x_vec,
      se = se_sigma_x_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    sx <- NULL
  }
  if ("tau_sqr_vec" %in% mat_names) {
    tau_sqr_vec <- OpenMx::mxEvalByName(
      name = "tau_sqr_vec",
      model = object$output,
      compute = TRUE
    )
    names_tau_sqr_vec <- rownames(tau_sqr_vec)
    tau_sqr_vec <- c(tau_sqr_vec)
    names(tau_sqr_vec) <- names_tau_sqr_vec
    se_tau_sqr_vec <- c(
      OpenMx::mxSE(
        x = "tau_sqr_vec",
        model = object$output,
        silent = TRUE
      )
    )
    t2 <- .CIWald(
      est = tau_sqr_vec,
      se = se_tau_sqr_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    t2 <- NULL
  }
  if ("i_sqr_vec" %in% mat_names) {
    i_sqr_vec <- OpenMx::mxEvalByName(
      name = "i_sqr_vec",
      model = object$output,
      compute = TRUE
    )
    names_i_sqr_vec <- rownames(i_sqr_vec)
    i_sqr_vec <- c(i_sqr_vec)
    names(i_sqr_vec) <- names_i_sqr_vec
    se_i_sqr_vec <- c(
      OpenMx::mxSE(
        x = "i_sqr_vec",
        model = object$output,
        silent = TRUE
      )
    )
    i2 <- .CIWald(
      est = i_sqr_vec,
      se = se_i_sqr_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    i2 <- NULL
  }
  if ("direct_vec" %in% mat_names) {
    direct_vec <- OpenMx::mxEvalByName(
      name = "direct_vec",
      model = object$output,
      compute = TRUE
    )
    names_direct_vec <- rownames(direct_vec)
    direct_vec <- c(direct_vec)
    names(direct_vec) <- names_direct_vec
    se_direct_vec <- c(
      OpenMx::mxSE(
        x = "direct_vec",
        model = object$output,
        silent = TRUE
      )
    )
    direct_vec <- .CIWald(
      est = direct_vec,
      se = se_direct_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    direct_vec <- NULL
  }
  if ("indirect_vec" %in% mat_names) {
    indirect_vec <- OpenMx::mxEvalByName(
      name = "indirect_vec",
      model = object$output,
      compute = TRUE
    )
    names_indirect_vec <- rownames(indirect_vec)
    indirect_vec <- c(indirect_vec)
    names(indirect_vec) <- names_indirect_vec
    se_indirect_vec <- c(
      OpenMx::mxSE(
        x = "indirect_vec",
        model = object$output,
        silent = TRUE
      )
    )
    indirect_vec <- .CIWald(
      est = indirect_vec,
      se = se_indirect_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    indirect_vec <- NULL
  }
  if ("total_vec" %in% mat_names) {
    total_vec <- OpenMx::mxEvalByName(
      name = "total_vec",
      model = object$output,
      compute = TRUE
    )
    names_total_vec <- rownames(total_vec)
    total_vec <- c(total_vec)
    names(total_vec) <- names_total_vec
    se_total_vec <- c(
      OpenMx::mxSE(
        x = "total_vec",
        model = object$output,
        silent = TRUE
      )
    )
    total_vec <- .CIWald(
      est = total_vec,
      se = se_total_vec,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    total_vec <- NULL
  }
  ie_xyz <- grep(
    "^ie_x\\d+_y\\d+_z\\d+$",
    mat_names,
    value = TRUE
  )
  if (length(ie_xyz) > 0) {
    ie_xyz <- do.call(
      what = "rbind",
      args = lapply(
        X = ie_xyz,
        FUN = function(i) {
          ie_xyz_i <- OpenMx::mxEvalByName(
            name = i,
            model = object$output,
            compute = TRUE
          )
          ie_xyz_i <- c(ie_xyz_i)
          names(ie_xyz_i) <- i
          se_ie_xyz_i <- c(
            OpenMx::mxSE(
              x = i,
              model = object$output,
              silent = TRUE
            )
          )
          ie_xyz_i <- .CIWald(
            est = ie_xyz_i,
            se = se_ie_xyz_i,
            theta = 0,
            alpha = alpha,
            z = TRUE,
            test = FALSE
          )
        }
      )
    )
  } else {
    ie_xyz <- NULL
  }
  ci <- list(
    y0 = y0,
    y1 = y1,
    yy = yy,
    z0 = z0,
    z1 = z1,
    zx = zx,
    psi = psi,
    x0 = x0,
    sx = sx,
    t2 = t2,
    i2 = i2,
    direct_vec = direct_vec,
    indirect_vec = indirect_vec,
    total_vec = total_vec,
    ie_xyz = ie_xyz
  )
  ci <- ci[
    !sapply(
      X = ci,
      FUN = is.null
    )
  ]
  do.call(
    what = "rbind",
    args = ci
  )
}
