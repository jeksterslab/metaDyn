.MetaTauSqr <- function(p,
                        tau_sqr_diag,
                        tau_sqr_d_free,
                        tau_sqr_d_values,
                        tau_sqr_d_lbound,
                        tau_sqr_d_ubound,
                        tau_sqr_l_free,
                        tau_sqr_l_values,
                        tau_sqr_l_lbound,
                        tau_sqr_l_ubound) {
  if (is.null(tau_sqr_d_lbound)) {
    tau_sqr_d_lbound <- -30
  }
  if (is.null(tau_sqr_d_ubound)) {
    tau_sqr_d_ubound <- 650
  }
  idx <- seq_len(p)
  if (tau_sqr_diag) {
    out <- .MetaTauSqrDiag(
      p = p,
      idx = idx,
      tau_sqr_d_free = tau_sqr_d_free,
      tau_sqr_d_values = tau_sqr_d_values,
      tau_sqr_d_lbound = tau_sqr_d_lbound,
      tau_sqr_d_ubound = tau_sqr_d_ubound,
      tau_sqr_d_equal = FALSE
    )
  } else {
    out <- .MetaTauSqrSym(
      p = p,
      idx = idx,
      tau_sqr_d_free = tau_sqr_d_free,
      tau_sqr_d_values = tau_sqr_d_values,
      tau_sqr_d_lbound = tau_sqr_d_lbound,
      tau_sqr_d_ubound = tau_sqr_d_ubound,
      tau_sqr_l_free = tau_sqr_l_free,
      tau_sqr_l_values = tau_sqr_l_values,
      tau_sqr_l_lbound = tau_sqr_l_lbound,
      tau_sqr_l_ubound = tau_sqr_l_ubound,
      tau_sqr_d_equal = FALSE
    )
  }
  out
}
