.MetaPsi <- function(r,
                     psi_diag,
                     psi_d_free,
                     psi_d_values,
                     psi_d_lbound,
                     psi_d_ubound,
                     psi_l_free,
                     psi_l_values,
                     psi_l_lbound,
                     psi_l_ubound) {
  if (is.null(psi_d_lbound)) {
    psi_d_lbound <- -30
  }
  if (is.null(psi_d_ubound)) {
    psi_d_ubound <- 650
  }
  idx <- seq_len(r)
  if (psi_diag) {
    out <- .MetaPsiDiag(
      r = r,
      idx = idx,
      psi_d_free = psi_d_free,
      psi_d_values = psi_d_values,
      psi_d_lbound = psi_d_lbound,
      psi_d_ubound = psi_d_ubound,
      psi_d_equal = FALSE
    )
  } else {
    out <- .MetaPsiSym(
      r = r,
      idx = idx,
      psi_d_free = psi_d_free,
      psi_d_values = psi_d_values,
      psi_d_lbound = psi_d_lbound,
      psi_d_ubound = psi_d_ubound,
      psi_l_free = psi_l_free,
      psi_l_values = psi_l_values,
      psi_l_lbound = psi_l_lbound,
      psi_l_ubound = psi_l_ubound,
      psi_d_equal = FALSE
    )
  }
  out
}
