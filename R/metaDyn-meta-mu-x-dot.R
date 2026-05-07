.MetaMuX <- function(m,
                     xnames,
                     mu_x_free,
                     mu_x_values,
                     mu_x_lbound,
                     mu_x_ubound) {
  .MxHelperFullMxMatrix(
    m = m,
    n = 1,
    free_val = mu_x_free,
    values = mu_x_values,
    lbound_val = mu_x_lbound,
    ubound_val = mu_x_ubound,
    vec = TRUE,
    row = xnames,
    col = "mu_x",
    name = "mu_x"
  )
}
