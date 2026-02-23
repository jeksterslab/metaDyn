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
  close_to <- function(val,
                       bnd) {
    na_bnd <- is.na(bnd) | !is.finite(bnd)
    na_val <- is.na(val) | !is.finite(val)
    sc <- pmax(
      abs(val),
      abs(bnd),
      1
    )
    tol <- abs_bnd_tol + rel_bnd_tol * sc
    res <- (
      !na_bnd
    ) & (
      !na_val
    ) & (
      abs(
        val - bnd
      ) <= tol
    )
    res[is.na(res)] <- FALSE
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
    any = any(
      atlb | atub,
      na.rm = TRUE
    ),
    at_lb = atlb,
    at_ub = atub
  )
}
