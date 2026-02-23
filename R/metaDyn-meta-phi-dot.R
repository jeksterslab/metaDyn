.MetaPhi <- function(p,
                     r,
                     ynames,
                     znames,
                     phi_free,
                     phi_values,
                     phi_lbound,
                     phi_ubound,
                     alpha) {
  .MxHelperFullMxMatrix(
    m = r,
    n = p,
    free_val = phi_free,
    values = phi_values,
    lbound_val = phi_lbound,
    ubound_val = phi_ubound,
    vec = TRUE,
    row = znames,
    col = ynames,
    name = "phi"
  )
}
