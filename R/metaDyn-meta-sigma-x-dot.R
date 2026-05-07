.MetaSigmaX <- function(m,
                        sigma_x_d_free,
                        sigma_x_d_values,
                        sigma_x_d_lbound,
                        sigma_x_d_ubound,
                        sigma_x_l_free,
                        sigma_x_l_values,
                        sigma_x_l_lbound,
                        sigma_x_l_ubound) {
  name <- "sigma_x"
  sigma_x_d <- paste0(
    name,
    "_d"
  )
  sigma_x_l <- paste0(
    name,
    "_l"
  )
  idx <- seq_len(m)
  .MxHelperSigmaFromLDLMxMatrix(
    p = m,
    name = name,
    column_name = sigma_x_d,
    sdiag_name = sigma_x_l,
    iden_name = NULL,
    d_free = sigma_x_d_free,
    d_values = sigma_x_d_values,
    d_lbound = sigma_x_d_lbound,
    d_ubound = sigma_x_d_ubound,
    d_rows = idx,
    d_cols = sigma_x_d,
    d_equal = FALSE,
    l_free = sigma_x_l_free,
    l_values = sigma_x_l_values,
    l_lbound = sigma_x_l_lbound,
    l_ubound = sigma_x_l_ubound,
    l_rows = idx,
    l_cols = idx
  )
}
