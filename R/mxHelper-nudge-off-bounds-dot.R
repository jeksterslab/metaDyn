.MxHelperNudgeOffBounds <- function(x,
                                    eps = 1e-6) {
  pars <- OpenMx::omxGetParameters(x)
  lb <- OpenMx::omxGetParameters(
    model = x,
    indep = FALSE,
    free = TRUE,
    fetch = "lbound"
  )
  ub <- OpenMx::omxGetParameters(
    model = x,
    indep = FALSE,
    free = TRUE,
    fetch = "ubound"
  )
  lb <- lb[names(pars)]
  ub <- ub[names(pars)]
  at_lb <- (
    !is.na(lb)
  ) & (
    is.finite(lb)
  ) & (
    abs(pars - lb) < eps
  )
  at_ub <- (
    !is.na(ub)
  ) & (
    is.finite(ub)
  ) & (
    abs(ub - pars) < eps
  )
  pars[at_lb] <- pars[at_lb] + 10 * eps
  pars[at_ub] <- pars[at_ub] - 10 * eps
  OpenMx::omxSetParameters(
    model = x,
    names(pars),
    values = pars
  )
}
