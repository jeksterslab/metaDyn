.MxHelperSdiagLabels <- function(p,
                                 name,
                                 sep) {
  out <- .MxHelperFullLabels(
    m = p,
    n = p,
    name = name,
    sep = sep
  )
  out[
    upper.tri(
      x = out,
      diag = TRUE
    )
  ] <- NA
  out
}
