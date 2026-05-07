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
      .MetaExpectedCovariancesStochasticXDistal(
        etaeta = etaeta,
        xnames = xnames,
        ynames = ynames,
        znames = znames
      )
    } else {
      .MetaExpectedCovariancesStochasticXNoDistal(
        etaeta = etaeta,
        xnames = xnames,
        ynames = ynames
      )
    }
  } else {
    if (distal) {
      .MetaExpectedCovariancesFixedXDistal(
        etaeta = etaeta,
        ynames = ynames,
        znames = znames
      )
    } else {
      .MetaExpectedCovariancesFixedXNoDistal(
        etaeta = etaeta,
        ynames = ynames
      )
    }
  }
}
