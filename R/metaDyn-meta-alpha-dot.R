.MetaAlpha <- function(p,
                       ynames,
                       alpha_free,
                       alpha_values,
                       alpha_lbound,
                       alpha_ubound,
                       alpha) {
  c(
    .MxHelperFullMxMatrix(
      m = p,
      n = 1,
      free_val = alpha_free,
      values = alpha_values,
      lbound_val = alpha_lbound,
      ubound_val = alpha_ubound,
      vec = TRUE,
      row = ynames,
      col = "alpha",
      name = "alpha"
    ),
    list(
      OpenMx::mxCI(
        reference = "alpha",
        interval = 1 - alpha
      )
    )
  )
}
