.MetaKappa <- function(r,
                       znames,
                       kappa_free,
                       kappa_values,
                       kappa_lbound,
                       kappa_ubound,
                       alpha) {
  .MxHelperFullMxMatrix(
    m = r,
    n = 1,
    free_val = kappa_free,
    values = kappa_values,
    lbound_val = kappa_lbound,
    ubound_val = kappa_ubound,
    vec = TRUE,
    row = znames,
    col = "kappa",
    name = "kappa"
  )
}
