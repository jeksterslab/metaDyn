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
  raw_data <- .PrepData(
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
  variables <- .MetaMxMatricesVariables(
    x = x,
    v = v,
    p = p,
    m = m,
    n = n,
    xnames = xnames,
    vnames = vnames,
    random = random,
    covariate = covariate
  )
  starts <- .MetaStarts(
    raw_data = raw_data,
    xnames = xnames,
    ynames = ynames,
    znames = znames,
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
  alpha_values <- starts$alpha_values
  tau_sqr_d_values <- starts$tau_sqr_d_values
  tau_sqr_l_values <- starts$tau_sqr_l_values
  gamma_values <- starts$gamma_values
  kappa_values <- starts$kappa_values
  phi_values <- starts$phi_values
  omega_values <- starts$omega_values
  psi_d_values <- starts$psi_d_values
  psi_l_values <- starts$psi_l_values
  parameters <- .MetaMxMatricesParameters(
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
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    gamma_lbound = gamma_lbound,
    gamma_ubound = gamma_ubound,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    kappa_lbound = kappa_lbound,
    kappa_ubound = kappa_ubound,
    phi_free = phi_free,
    phi_values = phi_values,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_free = omega_free,
    omega_values = omega_values,
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
  OpenMx::mxModel(
    model = "Model",
    variables,
    parameters,
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
      observed = raw_data
    ),
    OpenMx::mxExpectationNormal(
      covariance = "expected_covariance",
      means = "expected_mean",
      dimnames = dimnames
    ),
    OpenMx::mxFitFunctionML()
  )
}
