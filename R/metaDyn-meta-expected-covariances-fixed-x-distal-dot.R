.MetaExpectedCovariancesFixedXDistal <- function(etaeta,
                                                 ynames,
                                                 znames) {
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
      name = "expected_covariance_yy",
      dimnames = list(
        ynames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "phi %*% expected_covariance_etaeta",
      name = "expected_covariance_zy",
      dimnames = list(
        znames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "expected_covariance_etaeta %*% t(phi)",
      name = "expected_covariance_yz",
      dimnames = list(
        ynames,
        znames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "phi %*% expected_covariance_etaeta %*% t(phi)",
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
}
