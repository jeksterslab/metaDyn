.VHat <- function(v,
                  p,
                  n,
                  univariate) {
  if (univariate) {
    out <- .VHatUnivariate(
      v = v,
      p = p,
      n = n
    )
  } else {
    out <- .VHatMultivariate(
      v = v,
      p = p,
      n = n
    )
  }
  out
}
