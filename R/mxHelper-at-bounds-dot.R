.MxHelperAtBounds <- function(x,
                              abs_bnd_tol = 1e-6,
                              rel_bnd_tol = 1e-4,
                              allowed_bounds = NULL) {
  pars <- OpenMx::omxGetParameters(x)
  nm <- names(pars)

  # Bounds aligned to pars by name.
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

  # Normalize and validate allowed bounds.
  if (is.null(allowed_bounds)) {
    allowed_bounds <- list(
      lower = character(0),
      upper = character(0)
    )
  } else {
    if (!is.list(allowed_bounds)) {
      stop(
        "allowed_bounds must be NULL or a list."
      )
    }

    allowed_names <- names(allowed_bounds)

    if (
      !is.null(allowed_names) &&
        any(
          !allowed_names %in% c(
            "lower",
            "upper"
          )
        )
    ) {
      stop(
        paste0(
          "allowed_bounds may only contain ",
          "'lower' and 'upper'."
        )
      )
    }

    if (is.null(allowed_bounds$lower)) {
      allowed_bounds$lower <- character(0)
    }

    if (is.null(allowed_bounds$upper)) {
      allowed_bounds$upper <- character(0)
    }

    if (
      !is.character(allowed_bounds$lower) ||
        anyNA(allowed_bounds$lower)
    ) {
      stop(
        paste0(
          "allowed_bounds$lower must be ",
          "a character vector without missing values."
        )
      )
    }

    if (
      !is.character(allowed_bounds$upper) ||
        anyNA(allowed_bounds$upper)
    ) {
      stop(
        paste0(
          "allowed_bounds$upper must be ",
          "a character vector without missing values."
        )
      )
    }

    allowed_bounds$lower <- unique(
      allowed_bounds$lower
    )

    allowed_bounds$upper <- unique(
      allowed_bounds$upper
    )

    unknown <- setdiff(
      unique(
        c(
          allowed_bounds$lower,
          allowed_bounds$upper
        )
      ),
      nm
    )

    if (length(unknown) > 0L) {
      stop(
        paste0(
          "Unknown parameter label(s) in allowed_bounds: ",
          paste(
            unknown,
            collapse = ", "
          ),
          "."
        )
      )
    }
  }

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

  # Raw bound status is retained.
  actionable_lb <- atlb
  actionable_ub <- atub

  # Allowed bound hits remain visible in at_lb / at_ub,
  # but they are not actionable.
  if (length(allowed_bounds$lower) > 0L) {
    actionable_lb[
      nm %in% allowed_bounds$lower
    ] <- FALSE
  }

  if (length(allowed_bounds$upper) > 0L) {
    actionable_ub[
      nm %in% allowed_bounds$upper
    ] <- FALSE
  }

  list(
    any = any(
      atlb | atub,
      na.rm = TRUE
    ),
    actionable = any(
      actionable_lb | actionable_ub,
      na.rm = TRUE
    ),
    at_lb = atlb,
    at_ub = atub,
    actionable_lb = actionable_lb,
    actionable_ub = actionable_ub
  )
}
