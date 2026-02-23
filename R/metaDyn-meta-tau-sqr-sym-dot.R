.MetaTauSqrSym <- function(p,
                           idx,
                           tau_sqr_d_free,
                           tau_sqr_d_values,
                           tau_sqr_d_lbound,
                           tau_sqr_d_ubound,
                           tau_sqr_d_equal,
                           tau_sqr_l_fixed,
                           tau_sqr_l_free,
                           tau_sqr_l_values,
                           tau_sqr_l_lbound,
                           tau_sqr_l_ubound,
                           alpha) {
  name <- "tau_sqr"
  tau_sqr_d <- paste0(
    name,
    "_d"
  )
  tau_sqr_l <- paste0(
    name,
    "_l"
  )
  .MxHelperSigmaFromLDLMxMatrix(
    p = p,
    name = name,
    column_name = tau_sqr_d,
    sdiag_name = tau_sqr_l,
    iden_name = NULL,
    d_free = tau_sqr_d_free,
    d_values = tau_sqr_d_values,
    d_lbound = tau_sqr_d_lbound,
    d_ubound = tau_sqr_d_ubound,
    d_rows = idx,
    d_cols = tau_sqr_d,
    d_equal = tau_sqr_d_equal,
    l_free = tau_sqr_l_free,
    l_values = tau_sqr_l_values,
    l_lbound = tau_sqr_l_lbound,
    l_ubound = tau_sqr_l_ubound,
    l_rows = idx,
    l_cols = idx
  )
}
