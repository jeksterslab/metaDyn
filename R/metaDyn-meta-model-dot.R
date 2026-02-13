.MetaModel <- function(y,
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
                       alpha) {
  if (distal) {
    znames <- paste0(
      "z",
      seq_len(r)
    )
    dimnames <- znames
  } else {
    znames <- NULL
    dimnames <- c()
  }
  v_hat_univariate <- i_sqr_univariate
  idx <- seq_len(p)
  ynames <- paste0(
    "y",
    idx
  )
  dimnames <- c(
    dimnames,
    ynames
  )
  vnames <- .Vech(
    x = outer(
      X = idx,
      Y = idx,
      FUN = function(x, y) {
        paste0(
          "v",
          x,
          y
        )
      }
    )
  )
  if (covariate) {
    xnames <- paste0(
      "x",
      seq_len(m)
    )
  } else {
    xnames <- NULL
  }
  OpenMx::mxModel(
    model = "Model",
    .MetaMxMatricesVariables(
      x = x,
      v = v,
      p = p,
      m = m,
      n = n,
      xnames = xnames,
      vnames = vnames,
      random = random,
      covariate = covariate,
      v_hat_univariate = v_hat_univariate
    ),
    .MetaMxMatricesParameters(
      p = p,
      m = m,
      r = r,
      ynames = ynames,
      xnames = xnames,
      znames = znames,
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
      psi_l_ubound = psi_l_ubound,
      alpha = alpha
    ),
    .MetaExpectedMeans(
      covariate = covariate,
      distal = distal,
      ynames = ynames,
      znames = znames
    ),
    .MetaExpectedCovariances(
      random = random,
      distal = distal,
      ynames = ynames,
      znames = znames
    ),
    OpenMx::mxData(
      type = "raw",
      observed = .PrepData(
        y = y,
        v = v,
        x = x,
        z = z,
        ynames = ynames,
        vnames = vnames,
        xnames = xnames,
        znames = znames,
        covariate = covariate,
        distal = distal
      )
    ),
    OpenMx::mxExpectationNormal(
      covariance = "expected_covariance",
      means = "expected_mean",
      dimnames = dimnames
    ),
    OpenMx::mxFitFunctionML()
  )
}
