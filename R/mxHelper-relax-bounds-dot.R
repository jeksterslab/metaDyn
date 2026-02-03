.MxHelperRelaxBounds <- function(x,
                                 factor = 10) {
  lb <- OpenMx::omxGetParameters(
    model = x,
    fetch = "lbound"
  )
  ub <- OpenMx::omxGetParameters(
    model = x,
    fetch = "ubound"
  )
  nm <- names(lb)
  new_lb <- lb
  new_ub <- ub
  for (k in nm) {
    if (is.finite(new_lb[k])) {
      new_lb[k] <- new_lb[k] / factor
    }
    if (is.finite(new_ub[k])) {
      new_ub[k] <- new_ub[k] * factor
    }
  }
  OpenMx::omxSetParameters(
    model = x,
    labels = nm,
    lbound = new_lb,
    ubound = new_ub
  )
}
