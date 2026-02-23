.MetaTauSqrDiag <- function(p,
                            idx,
                            tau_sqr_d_free,
                            tau_sqr_d_values,
                            tau_sqr_d_lbound,
                            tau_sqr_d_ubound,
                            tau_sqr_d_equal,
                            alpha) {
  name <- "tau_sqr"
  tau_sqr_d <- paste0(
    name,
    "_d"
  )
  .MxHelperSigmaDiagFromLDLMxMatrix(
    p = p,
    name = name,
    column_name = tau_sqr_d,
    d_free = tau_sqr_d_free,
    d_values = tau_sqr_d_values,
    d_lbound = tau_sqr_d_lbound,
    d_ubound = tau_sqr_d_ubound,
    d_rows = idx,
    d_cols = tau_sqr_d,
    d_equal = tau_sqr_d_equal
  )
}
