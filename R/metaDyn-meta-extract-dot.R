.MetaExtract <- function(object) {
  p <- object$args$p
  m <- object$args$m
  r <- object$args$r
  coefs <- summary(object)
  parnames <- rownames(
    coefs
  )
  coefs <- coefs[, 1]
  names(coefs) <- parnames
  out <- list()
  alpha_names <- parnames[
    grep(
      pattern = "^alpha",
      x = parnames
    )
  ]
  if (length(alpha_names) > 0) {
    out <- c(
      out,
      alpha = list(unname(coefs[alpha_names]))
    )
  }
  beta_names <- parnames[
    grep(
      pattern = "^beta",
      x = parnames
    )
  ]
  if (length(beta_names) > 0) {
    beta <- matrix(
      data = 0,
      nrow = p,
      ncol = m
    )
    for (j in seq_len(m)) {
      for (i in seq_len(p)) {
        beta_ij <- paste0("beta[", i, ",", j, "]")
        if (
          paste0("beta[", i, ",", j, "]") %in% beta_names
        ) {
          beta[i, j] <- coefs[beta_ij]
        }
      }
    }
    out <- c(
      out,
      beta = list(beta)
    )
  }
  gamma_names <- parnames[
    grep(
      pattern = "^gamma",
      x = parnames
    )
  ]
  if (length(gamma_names) > 0) {
    gamma <- matrix(
      data = 0,
      nrow = p,
      ncol = m
    )
    for (j in seq_len(m)) {
      for (i in seq_len(p)) {
        gamma_ij <- paste0("gamma[", i, ",", j, "]")
        if (
          paste0("gamma[", i, ",", j, "]") %in% gamma_names
        ) {
          gamma[i, j] <- coefs[gamma_ij]
        }
      }
    }
    out <- c(
      out,
      gamma = list(gamma)
    )
  }
  tau_sqr_names <- parnames[
    grep(
      pattern = "^tau_sqr",
      x = parnames
    )
  ]
  if (length(tau_sqr_names) > 0) {
    tau_sqr <- matrix(
      data = 0,
      nrow = p,
      ncol = p
    )
    for (j in seq_len(p)) {
      for (i in seq_len(p)) {
        tau_sqr_ij <- paste0("tau_sqr[", i, ",", j, "]")
        if (
          paste0("tau_sqr[", i, ",", j, "]") %in% tau_sqr_names
        ) {
          tau_sqr[i, j] <- coefs[tau_sqr_ij]
          tau_sqr[j, i] <- tau_sqr[i, j]
        }
      }
    }
    out <- c(
      out,
      tau_sqr = list(tau_sqr)
    )
  }
  i_sqr_names <- parnames[
    grep(
      pattern = "^i_sqr",
      x = parnames
    )
  ]
  if (length(i_sqr_names) > 0) {
    out <- c(
      out,
      i_sqr = list(unname(coefs[i_sqr_names]))
    )
  }
  kappa_names <- parnames[
    grep(
      pattern = "^kappa",
      x = parnames
    )
  ]
  if (length(kappa_names) > 0) {
    out <- c(
      out,
      kappa = list(unname(coefs[kappa_names]))
    )
  }
  phi_names <- parnames[
    grep(
      pattern = "^phi",
      x = parnames
    )
  ]
  if (length(phi_names) > 0) {
    phi <- matrix(
      data = 0,
      nrow = r,
      ncol = p
    )
    for (j in seq_len(p)) {
      for (i in seq_len(r)) {
        phi_ij <- paste0("phi[", i, ",", j, "]")
        if (
          paste0("phi[", i, ",", j, "]") %in% phi_names
        ) {
          phi[i, j] <- coefs[phi_ij]
        }
      }
    }
    out <- c(
      out,
      phi = list(phi)
    )
  }
  psi_names <- parnames[
    grep(
      pattern = "^psi",
      x = parnames
    )
  ]
  if (length(psi_names) > 0) {
    psi <- matrix(
      data = 0,
      nrow = r,
      ncol = r
    )
    for (j in seq_len(r)) {
      for (i in seq_len(r)) {
        psi_ij <- paste0("psi[", i, ",", j, "]")
        if (
          paste0("psi[", i, ",", j, "]") %in% psi_names
        ) {
          psi[i, j] <- coefs[psi_ij]
          psi[j, i] <- psi[i, j]
        }
      }
    }
    out <- c(
      out,
      psi = list(psi)
    )
  }
  omega_names <- parnames[
    grep(
      pattern = "^omega",
      x = parnames
    )
  ]
  if (length(omega_names) > 0) {
    omega <- matrix(
      data = 0,
      nrow = r,
      ncol = m
    )
    for (j in seq_len(m)) {
      for (i in seq_len(r)) {
        omega_ij <- paste0("omega[", i, ",", j, "]")
        if (
          paste0("omega[", i, ",", j, "]") %in% omega_names
        ) {
          omega[i, j] <- coefs[omega_ij]
        }
      }
    }
    out <- c(
      out,
      omega = list(omega)
    )
  }
  out
}
