.MetaGamma <- function(p,
                       m,
                       ynames,
                       xnames,
                       gamma_free,
                       gamma_values,
                       gamma_lbound,
                       gamma_ubound,
                       alpha) {
  .MxHelperFullMxMatrix(
    m = p,
    n = m,
    free_val = gamma_free,
    values = gamma_values,
    lbound_val = gamma_lbound,
    ubound_val = gamma_ubound,
    vec = TRUE,
    row = ynames,
    col = xnames,
    name = "gamma"
  )
}
