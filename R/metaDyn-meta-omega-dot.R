.MetaOmega <- function(m,
                       r,
                       xnames,
                       znames,
                       omega_free,
                       omega_values,
                       omega_lbound,
                       omega_ubound) {
  .MxHelperFullMxMatrix(
    m = r,
    n = m,
    free_val = omega_free,
    values = omega_values,
    lbound_val = omega_lbound,
    ubound_val = omega_ubound,
    vec = TRUE,
    row = znames,
    col = xnames,
    name = "omega"
  )
}
