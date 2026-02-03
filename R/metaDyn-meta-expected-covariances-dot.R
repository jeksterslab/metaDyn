.MetaExpectedCovariances <- function(random,
                                     distal,
                                     ynames,
                                     znames) {
  if (random) {
    yy <- "v + tau_sqr"
  } else {
    yy <- "v"
  }
  if (distal) {
    out <- list(
      OpenMx::mxAlgebraFromString(
        algString = yy,
        name = "expected_covariance_yy",
        dimnames = list(
          ynames,
          ynames
        )
      ),
      OpenMx::mxAlgebraFromString(
        algString = "phi %*% expected_covariance_yy",
        name = "expected_covariance_zy",
        dimnames = list(
          znames,
          ynames
        )
      ),
      OpenMx::mxAlgebraFromString(
        algString = "expected_covariance_yy %*% t(phi)",
        name = "expected_covariance_yz",
        dimnames = list(
          ynames,
          znames
        )
      ),
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "phi %*% expected_covariance_yy %*% t(phi)",
          " + ",
          "psi"
        ),
        name = "expected_covariance_zz",
        dimnames = list(
          znames,
          znames
        )
      ),
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "rbind(",
          "cbind(expected_covariance_zz, expected_covariance_zy)",
          ",",
          "cbind(expected_covariance_yz, expected_covariance_yy)",
          ")"
        ),
        name = "expected_covariance",
        dimnames = list(
          c(znames, ynames),
          c(znames, ynames)
        )
      )
    )
  } else {
    out <- list(
      OpenMx::mxAlgebraFromString(
        algString = yy,
        name = "expected_covariance",
        dimnames = list(
          ynames,
          ynames
        )
      )
    )
  }
  out
}
