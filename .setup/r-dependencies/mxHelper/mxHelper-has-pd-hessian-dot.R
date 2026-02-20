.MxHelperHasPdHessian <- function(x,
                                  hess_tol_abs = 1e-8,
                                  hess_tol_rel = 1e-10,
                                  check_condition = FALSE,
                                  cond_max = 1e12) {
  hessian <- x$output$hessian
  if (is.null(hessian) || !is.matrix(hessian) || anyNA(hessian)) {
    # nolint start
    return(FALSE)
    # nolint end
  }
  hessian <- 0.5 * (hessian + t(hessian))
  eig <- tryCatch(
    eigen(
      hessian,
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    error = function(e) {
      NA_real_
    }
  )
  if (!is.numeric(eig) || length(eig) == 0L || !all(is.finite(eig))) {
    # nolint start
    return(FALSE)
    # nolint end
  }
  scale <- max(1, max(abs(eig)))
  thresh <- hess_tol_abs + hess_tol_rel * scale
  is_pd <- min(eig) > thresh
  if (!is_pd) {
    # nolint start
    return(FALSE)
    # nolint end
  }
  if (check_condition) {
    # Condition number approx: max/min eigenvalue for symmetric PD matrix
    cond <- max(eig) / min(eig)
    if (!is.finite(cond) || cond > cond_max) {
      # nolint start
      return(FALSE)
      # nolint end
    }
  }
  TRUE
}
