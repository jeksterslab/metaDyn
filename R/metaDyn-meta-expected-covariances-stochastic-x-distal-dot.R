.MetaExpectedCovariancesStochasticXDistal <- function(etaeta,
                                                      xnames,
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
      algString = "phi %*% gamma + omega",
      name = "expected_covariance_zx_slope",
      dimnames = list(
        znames,
        xnames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "expected_covariance_zx_slope %*% expected_covariance_xx",
      name = "expected_covariance_zx",
      dimnames = list(
        znames,
        xnames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = "expected_covariance_xx %*% t(expected_covariance_zx_slope)",
      name = "expected_covariance_xz",
      dimnames = list(
        xnames,
        znames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "expected_covariance_zx_slope",
        " %*% expected_covariance_xx %*% t(gamma)",
        " + ",
        "phi %*% expected_covariance_etaeta"
      ),
      name = "expected_covariance_zy",
      dimnames = list(
        znames,
        ynames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "gamma %*% expected_covariance_xx",
        " %*% t(expected_covariance_zx_slope)",
        " + ",
        "expected_covariance_etaeta %*% t(phi)"
      ),
      name = "expected_covariance_yz",
      dimnames = list(
        ynames,
        znames
      )
    ),
    OpenMx::mxAlgebraFromString(
      algString = paste0(
        "expected_covariance_zx_slope",
        " %*% expected_covariance_xx",
        " %*% t(expected_covariance_zx_slope)",
        " + ",
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
        "cbind(expected_covariance_zz, expected_covariance_zy, expected_covariance_zx)",
        ",",
        "cbind(expected_covariance_yz, expected_covariance_yy, expected_covariance_yx)",
        ",",
        "cbind(expected_covariance_xz, expected_covariance_xy, expected_covariance_xx)",
        ")"
      ),
      name = "expected_covariance",
      dimnames = list(
        c(znames, ynames, xnames),
        c(znames, ynames, xnames)
      )
    )
  )
}
