.MxHelperRun <- function(model,
                         grad_tol = 1e-2,
                         ok_codes = 0L,
                         require_finite_fit = TRUE,
                         hess_tol_abs = 1e-8,
                         hess_tol_rel = 1e-10,
                         check_condition = FALSE,
                         cond_max = 1e12,
                         silent = FALSE,
                         ...) {
  obj <- inherits(model, "MxModel")
  if (isFALSE(obj)) {
    # nolint start
    return(NULL)
    # nolint end
  }
  # Default: rerun unless we can prove it's already "good"
  run <- TRUE
  has_status <- (
    !is.null(model$output) &&
      !is.null(model$output$status) &&
      !is.null(model$output$status$code)
  )
  if (has_status && model$output$status$code == 0L) {
    good_fit <- tryCatch(
      .MxHelperIsGoodFit(
        x = model,
        grad_tol = grad_tol,
        ok_codes = ok_codes,
        require_finite_fit = require_finite_fit
      ),
      error = function(e) FALSE
    )
    if (isTRUE(good_fit)) {
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
      run <- isFALSE(pd_hessian)
    } else {
      run <- TRUE
    }
  }
  input_model <- model
  if (run) {
    # allow status codes 0, 1, 5, 6
    model <- tryCatch(
      OpenMx::mxTryHard(
        model = model,
        silent = silent,
        OKstatuscodes = c(0, 1, 5, 6),
        ...
      ),
      error = function(e) {
        input_model
      }
    )
  }
  model
}
