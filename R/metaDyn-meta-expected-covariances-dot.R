.MetaExpectedCovariances <- function(random,
                                     distal,
                                     fixed_x,
                                     xnames,
                                     ynames,
                                     znames) {
  stochastic_x <- (
    !fixed_x &&
      !is.null(xnames)
  )
  if (random) {
    etaeta <- "tau_sqr"
  } else {
    etaeta <- "0 * v"
  }
  if (stochastic_x) {
    if (distal) {
      .MetaCovJointXZ(
        etaeta = etaeta,
        xnames = xnames,
        ynames = ynames,
        znames = znames
      )
    } else {
      .MetaCovJointXY(
        etaeta = etaeta,
        xnames = xnames,
        ynames = ynames
      )
    }
  } else {
    if (distal) {
      .MetaCovCondZ(
        etaeta = etaeta,
        ynames = ynames,
        znames = znames
      )
    } else {
      .MetaCovCondY(
        etaeta = etaeta,
        ynames = ynames
      )
    }
  }
}
