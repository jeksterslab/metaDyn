.MxHelperIsGoodFit <- function(x,
                               grad_tol = 1e-2,
                               ok_codes = 0L,
                               require_finite_fit = TRUE) {
  status_code <- tryCatch(
    {
      x$output$status$code
    },
    error = function(e) {
      NA_integer_
    }
  )
  if (!isTRUE(status_code %in% ok_codes)) {
    # nolint start
    return(FALSE)
    # nolint end
  }
  if (require_finite_fit) {
    fitval <- tryCatch(
      {
        x$output$fit
      },
      error = function(e) {
        NA_real_
      }
    )
    if (!is.finite(fitval)) {
      # nolint start
      return(FALSE)
      # nolint end
    }
  }
  grad_raw <- tryCatch(
    {
      x$output$gradient
    },
    error = function(e) {
      NULL
    }
  )
  if (is.numeric(grad_raw) && length(grad_raw) > 0) {
    max_grad <- suppressWarnings(max(abs(grad_raw), na.rm = TRUE))
    if (!is.finite(max_grad)) {
      # nolint start
      return(FALSE)
      # nolint end
    }
  } else {
    # Avoid summary unless needed
    max_grad <- tryCatch(
      {
        as.numeric(summary(x)$maxAbsGradient)
      },
      error = function(e) {
        NA_real_
      }
    )
  }
  is.finite(max_grad) && (max_grad < grad_tol)
}
