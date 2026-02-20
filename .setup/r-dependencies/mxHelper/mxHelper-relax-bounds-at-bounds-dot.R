.MxHelperRelaxBoundsAtBounds <- function(x,
                                         factor = 10,
                                         abs_bnd_tol = 1e-6,
                                         rel_bnd_tol = 1e-4,
                                         exclude = NULL,
                                         protect_lb_zero = TRUE) {
  if (
    !is.numeric(factor) ||
      length(factor) != 1L ||
      !is.finite(factor) ||
      factor <= 1
  ) {
    stop(
      ".MxHelperRelaxBoundsAtBounds() requires factor > 1."
    )
  }
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
  nm <- intersect(
    intersect(
      names(lb),
      names(ub)
    ),
    names(pars)
  )
  bd <- .MxHelperAtBounds(
    x = x,
    abs_bnd_tol = abs_bnd_tol,
    rel_bnd_tol = rel_bnd_tol
  )
  hits <- (
    bd$at_lb[nm] | bd$at_ub[nm]
  )
  hits[is.na(hits)] <- FALSE
  targets <- nm[hits]
  if (!is.null(exclude)) {
    targets <- setdiff(
      targets,
      exclude
    )
  }
  if (length(targets) == 0L) {
    # nolint start
    return(x)
    # nolint end
  }
  new_lb <- lb
  new_ub <- ub
  for (k in targets) {
    if (
      isTRUE(bd$at_lb[k]) && is.finite(lb[k])
    ) {
      if (
        !(protect_lb_zero && isTRUE(lb[k] == 0))
      ) {
        new_lb[k] <- if (lb[k] < 0) {
          lb[k] * factor
        } else {
          lb[k] / factor
        }
      }
    }

    if (isTRUE(bd$at_ub[k]) && is.finite(ub[k])) {
      new_ub[k] <- if (ub[k] > 0) {
        ub[k] * factor
      } else {
        ub[k] / factor
      }
    }

    if (
      is.finite(new_lb[k]) && is.finite(new_ub[k]) && new_lb[k] >= new_ub[k]
    ) {
      mid <- (
        new_lb[k] + new_ub[k]
      ) / 2
      span <- abs(mid) + 1
      new_lb[k] <- mid - span
      new_ub[k] <- mid + span
    }
  }
  OpenMx::omxSetParameters(
    model = x,
    labels = targets,
    lbound = new_lb[targets],
    ubound = new_ub[targets]
  )
}
