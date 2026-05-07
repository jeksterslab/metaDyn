.MetaCovCondY <- function(etaeta,
                          ynames) {
  list(
    OpenMx::mxAlgebraFromString(
      algString = etaeta,
      name = "expected_covariance_etaeta",
      dimnames = list(
        ynames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "expected_covariance_etaeta + v",
      name = "expected_covariance",
      dimnames = list(
        ynames,
        ynames
      )
    )
  )
}
