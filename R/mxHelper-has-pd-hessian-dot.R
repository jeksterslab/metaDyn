.MxHelperHasPdHessian <- function(x,
                                  tol = 1e-8) {
  hessian <- x$output$hessian
  if (is.null(hessian)) {
    out <- FALSE
  } else {
    eigen_values <- tryCatch(
      expr = eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    out <- (
      isTRUE(all(is.finite(eigen_values)))
    ) && (
      all(eigen_values > tol)
    )
  }
  out
}
