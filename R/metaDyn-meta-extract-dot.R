.MetaExtract <- function(object) {
  mat_names <- c(
    names(object$output$algebras),
    names(object$output$matrices)
  )
  if ("alpha" %in% mat_names) {
    alpha <- OpenMx::mxEvalByName(
      name = "alpha",
      model = object$output,
      compute = TRUE
    )
  } else {
    alpha <- NULL
  }
  if ("beta" %in% mat_names) {
    beta <- OpenMx::mxEvalByName(
      name = "beta",
      model = object$output,
      compute = TRUE
    )
  } else {
    beta <- NULL
  }
  if ("gamma" %in% mat_names) {
    gamma <- OpenMx::mxEvalByName(
      name = "gamma",
      model = object$output,
      compute = TRUE
    )
  } else {
    gamma <- NULL
  }
  if ("kappa" %in% mat_names) {
    kappa <- OpenMx::mxEvalByName(
      name = "kappa",
      model = object$output,
      compute = TRUE
    )
  } else {
    kappa <- NULL
  }
  if ("phi" %in% mat_names) {
    phi <- OpenMx::mxEvalByName(
      name = "phi",
      model = object$output,
      compute = TRUE
    )
  } else {
    phi <- NULL
  }
  if ("omega" %in% mat_names) {
    omega <- OpenMx::mxEvalByName(
      name = "omega",
      model = object$output,
      compute = TRUE
    )
  } else {
    omega <- NULL
  }
  if ("psi" %in% mat_names) {
    psi <- OpenMx::mxEvalByName(
      name = "psi",
      model = object$output,
      compute = TRUE
    )
    colnames(psi) <- rownames(psi) <- rownames(phi)
  } else {
    psi <- NULL
  }
  if ("tau_sqr" %in% mat_names) {
    tau_sqr <- OpenMx::mxEvalByName(
      name = "tau_sqr",
      model = object$output,
      compute = TRUE
    )
    colnames(tau_sqr) <- rownames(tau_sqr) <- rownames(alpha)
  } else {
    tau_sqr <- NULL
  }
  if ("i_sqr" %in% mat_names) {
    i_sqr <- OpenMx::mxEvalByName(
      name = "i_sqr",
      model = object$output,
      compute = TRUE
    )
    rownames(i_sqr) <- rownames(alpha)
  } else {
    i_sqr <- NULL
  }
  list(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    kappa = kappa,
    phi = phi,
    omega = omega,
    psi = psi,
    tau_sqr = tau_sqr,
    i_sqr = i_sqr
  )
}
