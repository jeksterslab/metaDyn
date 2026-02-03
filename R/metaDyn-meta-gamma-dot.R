.MetaGamma <- function(p,
                       m,
                       ynames,
                       xnames,
                       gamma_values,
                       gamma_free,
                       gamma_lbound,
                       gamma_ubound) {
  c(
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
    ),
    list(
      OpenMx::mxCI(
        reference = "gamma",
        interval = 0.95
      )
    )
  )
}
