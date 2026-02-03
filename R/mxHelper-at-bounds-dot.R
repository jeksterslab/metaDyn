.MxHelperAtBounds <- function(x,
                              abs_bnd_tol = 1e-6,
                              rel_bnd_tol = 1e-4) {
  pars <- OpenMx::omxGetParameters(x)
  nm <- names(pars)

  # bounds aligned to pars by name
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

  close_to <- function(val, bnd) {
    # NA bound -> never at that bound
    na <- is.na(bnd) | !is.finite(bnd)
    sc <- pmax(abs(val), abs(bnd), 1) # scale
    tol <- abs_bnd_tol + rel_bnd_tol * sc
    res <- !na & (abs(val - bnd) <= tol)
    res
  }

  atlb <- close_to(
    val = pars,
    bnd = lb
  )
  atub <- close_to(
    val = pars,
    bnd = ub
  )

  list(
    any = any(atlb | atub),
    at_lb = atlb,
    at_ub = atub
  )
}
