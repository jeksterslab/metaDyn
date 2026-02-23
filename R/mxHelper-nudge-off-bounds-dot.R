.MxHelperNudgeOffBounds <- function(x,
                                    abs_bnd_tol = 1e-6,
                                    rel_bnd_tol = 1e-4,
                                    mult = 10) {
  pars <- OpenMx::omxGetParameters(x)
  nm <- names(pars)
  lb <- OpenMx::omxGetParameters(
    model = x,
    indep = FALSE,
    free = TRUE,
    fetch = "lbound"
  )[nm]
  ub <- OpenMx::omxGetParameters(
    model = x,
    indep = FALSE,
    free = TRUE,
    fetch = "ubound"
  )[nm]
  bd <- .MxHelperAtBounds(
    x = x,
    abs_bnd_tol = abs_bnd_tol,
    rel_bnd_tol = rel_bnd_tol
  )
  # Scale used for step sizes; replace non-finite scale with 1
  sc_lb <- pmax(abs(pars), abs(lb), 1)
  sc_lb[!is.finite(sc_lb)] <- 1
  sc_ub <- pmax(abs(pars), abs(ub), 1)
  sc_ub[!is.finite(sc_ub)] <- 1
  # Step size that clears the tolerance band used by .MxHelperAtBounds()
  step_lb <- mult * (abs_bnd_tol + rel_bnd_tol * sc_lb)
  step_ub <- mult * (abs_bnd_tol + rel_bnd_tol * sc_ub)
  # Inner clamp amount (same style as AtBounds)
  inner_tol <- function(val, bnd) {
    sc <- pmax(abs(val), abs(bnd), 1)
    sc[!is.finite(sc)] <- 1
    abs_bnd_tol + rel_bnd_tol * sc
  }
  # Nudge inward from lower bounds
  if (any(bd$at_lb, na.rm = TRUE)) {
    idx <- which(bd$at_lb)
    pars[idx] <- pars[idx] + step_lb[idx]

    finite_ub <- is.finite(ub[idx])
    if (any(finite_ub)) {
      iu <- idx[finite_ub]
      pars[iu] <- pmin(
        pars[iu],
        ub[iu] - inner_tol(pars[iu], ub[iu])
      )
    }
  }
  # Nudge inward from upper bounds
  if (any(bd$at_ub, na.rm = TRUE)) {
    idx <- which(bd$at_ub)
    pars[idx] <- pars[idx] - step_ub[idx]

    finite_lb <- is.finite(lb[idx])
    if (any(finite_lb)) {
      il <- idx[finite_lb]
      pars[il] <- pmax(
        pars[il],
        lb[il] + inner_tol(pars[il], lb[il])
      )
    }
  }
  OpenMx::omxSetParameters(
    model = x,
    labels = nm,
    values = pars
  )
}
