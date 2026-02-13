.MetaPsiDiag <- function(r,
                         idx,
                         psi_d_free,
                         psi_d_values,
                         psi_d_lbound,
                         psi_d_ubound,
                         psi_d_equal,
                         alpha) {
  name <- "psi"
  psi_d <- paste0(
    name,
    "_d"
  )
  c(
    .MxHelperSigmaDiagFromLDLMxMatrix(
      p = r,
      name = name,
      column_name = psi_d,
      d_free = psi_d_free,
      d_values = psi_d_values,
      d_lbound = psi_d_lbound,
      d_ubound = psi_d_ubound,
      d_rows = idx,
      d_cols = psi_d,
      d_equal = psi_d_equal
    ),
    list(
      OpenMx::mxCI(
        reference = name,
        interval = 1 - alpha
      )
    )
  )
}
