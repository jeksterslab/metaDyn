.MxHelperIsGoodFit <- function(x,
                               tol = 1e-2) {
  status_code <- tryCatch(
    x$output$status$code,
    error = function(e) {
      NA_integer_
    }
  )
  grad_raw <- tryCatch(
    x$output$gradient,
    error = function(e) {
      NULL
    }
  )
  max_grad <- NA_real_
  if (is.numeric(grad_raw)) {
    max_grad <- suppressWarnings(max(abs(grad_raw)))
  } else {
    max_grad <- tryCatch(
      {
        s <- summary(x)
        as.numeric(s$maxAbsGradient)
      },
      error = function(e) {
        NA_real_
      }
    )
  }
  (
    isTRUE(status_code == 0L)
  ) && (
    is.finite(max_grad)
  ) && (
    max_grad < tol
  )
}
