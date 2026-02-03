.CIWaldMeta <- function(object,
                        alpha) {
  coef_alpha <- OpenMx::mxEvalByName(
    name = "alpha_vec",
    model = object$output,
    compute = TRUE
  )
  names_alpha <- rownames(coef_alpha)
  coef_alpha <- c(coef_alpha)
  names(coef_alpha) <- names_alpha
  se_alpha <- c(
    OpenMx::mxSE(
      x = "alpha_vec",
      model = object$output,
      silent = TRUE
    )
  )
  y0 <- .CIWald(
    est = coef_alpha,
    se = se_alpha,
    theta = 0,
    alpha = alpha,
    z = TRUE,
    test = FALSE
  )
  if (object$args$covariate) {
    coef_gamma <- OpenMx::mxEvalByName(
      name = "gamma_vec",
      model = object$output,
      compute = TRUE
    )
    names_gamma <- rownames(coef_gamma)
    coef_gamma <- c(coef_gamma)
    names(coef_gamma) <- names_gamma
    se_gamma <- c(
      OpenMx::mxSE(
        x = "gamma_vec",
        model = object$output,
        silent = TRUE
      )
    )
    y1 <- .CIWald(
      est = coef_gamma,
      se = se_gamma,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    y1 <- NULL
  }
  if (object$args$distal) {
    coef_kappa <- OpenMx::mxEvalByName(
      name = "kappa_vec",
      model = object$output,
      compute = TRUE
    )
    names_kappa <- rownames(coef_kappa)
    coef_kappa <- c(coef_kappa)
    names(coef_kappa) <- names_kappa
    se_kappa <- c(
      OpenMx::mxSE(
        x = "kappa_vec",
        model = object$output,
        silent = TRUE
      )
    )
    z0 <- .CIWald(
      est = coef_kappa,
      se = se_kappa,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
    coef_phi <- OpenMx::mxEvalByName(
      name = "phi_vec",
      model = object$output,
      compute = TRUE
    )
    names_phi <- rownames(coef_phi)
    coef_phi <- c(coef_phi)
    names(coef_phi) <- names_phi
    se_phi <- c(
      OpenMx::mxSE(
        x = "phi_vec",
        model = object$output,
        silent = TRUE
      )
    )
    z1 <- .CIWald(
      est = coef_phi,
      se = se_phi,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
    coef_psi <- OpenMx::mxEvalByName(
      name = "psi_vec",
      model = object$output,
      compute = TRUE
    )
    names_psi <- rownames(coef_psi)
    coef_psi <- c(coef_psi)
    names(coef_psi) <- names_psi
    se_psi <- c(
      OpenMx::mxSE(
        x = "psi_vec",
        model = object$output,
        silent = TRUE
      )
    )
    psi <- .CIWald(
      est = coef_psi,
      se = se_psi,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
    if (object$args$covariate) {
      coef_omega <- OpenMx::mxEvalByName(
        name = "omega_vec",
        model = object$output,
        compute = TRUE
      )
      names_omega <- rownames(coef_omega)
      coef_omega <- c(coef_omega)
      names(coef_omega) <- names_omega
      se_omega <- c(
        OpenMx::mxSE(
          x = "omega_vec",
          model = object$output,
          silent = TRUE
        )
      )
      zx <- .CIWald(
        est = coef_omega,
        se = se_omega,
        theta = 0,
        alpha = alpha,
        z = TRUE,
        test = FALSE
      )
    } else {
      zx <- NULL
    }
  } else {
    z0 <- NULL
    z1 <- NULL
    psi <- NULL
    zx <- NULL
  }
  if (object$args$random) {
    coef_tau_sqr <- OpenMx::mxEvalByName(
      name = "tau_sqr_vec",
      model = object$output,
      compute = TRUE
    )
    names_tau_sqr <- rownames(coef_tau_sqr)
    coef_tau_sqr <- c(coef_tau_sqr)
    names(coef_tau_sqr) <- names_tau_sqr
    se_tau_sqr <- c(
      OpenMx::mxSE(
        x = "tau_sqr_vec",
        model = object$output,
        silent = TRUE
      )
    )
    t2 <- .CIWald(
      est = coef_tau_sqr,
      se = se_tau_sqr,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
    coef_i_sqr <- OpenMx::mxEvalByName(
      name = "i_sqr_vec",
      model = object$output,
      compute = TRUE
    )
    names_i_sqr <- gsub(
      pattern = "^alpha",
      replacement = "i_sqr",
      x = names_alpha
    )
    coef_i_sqr <- c(coef_i_sqr)
    names(coef_i_sqr) <- names_i_sqr
    se_i_sqr <- c(
      OpenMx::mxSE(
        x = "i_sqr_vec",
        model = object$output,
        silent = TRUE
      )
    )
    i2 <- .CIWald(
      est = coef_i_sqr,
      se = se_i_sqr,
      theta = 0,
      alpha = alpha,
      z = TRUE,
      test = FALSE
    )
  } else {
    t2 <- NULL
    i2 <- NULL
  }
  ci <- list(
    y0 = y0,
    y1 = y1,
    z0 = z0,
    z1 = z1,
    zx = zx,
    psi = psi,
    t2 = t2,
    i2 = i2
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
