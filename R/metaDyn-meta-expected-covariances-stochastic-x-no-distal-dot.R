.MetaExpectedCovariancesStochasticXNoDistal <- function(etaeta,
                                                        xnames,
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
      algString = "sigma_x",
      name = "expected_covariance_xx",
      dimnames = list(
        xnames,
        xnames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "gamma %*% expected_covariance_xx %*% t(gamma)",
        " + ",
        "expected_covariance_etaeta",
        " + ",
        "v"
      ),
      name = "expected_covariance_yy",
      dimnames = list(
        ynames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "gamma %*% expected_covariance_xx",
      name = "expected_covariance_yx",
      dimnames = list(
        ynames,
        xnames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "expected_covariance_xx %*% t(gamma)",
      name = "expected_covariance_xy",
      dimnames = list(
        xnames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "rbind(",
        "cbind(expected_covariance_yy, expected_covariance_yx)",
        ",",
        "cbind(expected_covariance_xy, expected_covariance_xx)",
        ")"
      ),
      name = "expected_covariance",
      dimnames = list(
        c(ynames, xnames),
        c(ynames, xnames)
      )
    )
  )
}
