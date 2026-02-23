.MxHelperNeedsRescue <- function(model,
                                 grad_tol = 1e-2,
                                 ok_codes = 0L,
                                 require_finite_fit = TRUE,
                                 hess_tol_abs = 1e-8,
                                 hess_tol_rel = 1e-10,
                                 check_condition = FALSE,
                                 cond_max = 1e12,
                                 abs_bnd_tol = 1e-6,
                                 rel_bnd_tol = 1e-4) {
  bad_status <- (
    is.null(model$output) ||
      is.null(model$output$status) ||
      is.null(model$output$status$code) ||
      is.na(model$output$status$code) ||
      !(model$output$status$code %in% ok_codes)
  )

  if (isTRUE(bad_status)) {
    # nolint start
    return(TRUE)
    # nolint end
  }
  good_fit <- tryCatch(
    .MxHelperIsGoodFit(
      x = model,
      grad_tol = grad_tol,
      ok_codes = ok_codes,
      require_finite_fit = require_finite_fit
    ),
    error = function(e) FALSE
  )
  if (!isTRUE(good_fit)) {
    # nolint start
    return(TRUE)
    # nolint end
  }
  pd_hessian <- tryCatch(
    .MxHelperHasPdHessian(
      x = model,
      hess_tol_abs = hess_tol_abs,
      hess_tol_rel = hess_tol_rel,
      check_condition = check_condition,
      cond_max = cond_max
    ),
    error = function(e) FALSE
  )
  if (!isTRUE(pd_hessian)) {
    # nolint start
    return(TRUE)
    # nolint end
  }
  at_bounds_any <- tryCatch(
    .MxHelperAtBounds(
      x = model,
      abs_bnd_tol = abs_bnd_tol,
      rel_bnd_tol = rel_bnd_tol
    )$any,
    error = function(e) FALSE
  )
  isTRUE(at_bounds_any)
}
