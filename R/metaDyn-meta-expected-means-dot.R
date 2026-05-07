.MetaExpectedMeans <- function(covariate,
                               distal,
                               fixed_x,
                               xnames,
                               ynames,
                               znames) {
  if (covariate) {
    if (fixed_x) {
      expected_mean_y_alg <- "t(alpha + gamma %*% x)"
      expected_mean_z_alg <- paste0(
        "t(",
        "kappa",
        " + ",
        "phi %*% (alpha + gamma %*% x)",
        " + ",
        "(omega %*% x)",
        ")"
      )
    } else {
      expected_mean_y_alg <- "t(alpha + gamma %*% mu_x)"
      expected_mean_z_alg <- paste0(
        "t(",
        "kappa",
        " + ",
        "phi %*% (alpha + gamma %*% mu_x)",
        " + ",
        "(omega %*% mu_x)",
        ")"
      )
    }

    if (distal) {
      if (fixed_x) {
        out <- list(
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_z_alg,
            name = "expected_mean_z",
            dimnames = list(
              "1",
              znames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_y_alg,
            name = "expected_mean_y",
            dimnames = list(
              "1",
              ynames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = "cbind(expected_mean_z, expected_mean_y)",
            name = "expected_mean",
            dimnames = list(
              "1",
              c(znames, ynames)
            )
          )
        )
      } else {
        out <- list(
          OpenMx::mxAlgebraFromString(
            algString = "t(mu_x)",
            name = "expected_mean_x",
            dimnames = list(
              "1",
              xnames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_z_alg,
            name = "expected_mean_z",
            dimnames = list(
              "1",
              znames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_y_alg,
            name = "expected_mean_y",
            dimnames = list(
              "1",
              ynames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = "cbind(expected_mean_z, expected_mean_y, expected_mean_x)",
            name = "expected_mean",
            dimnames = list(
              "1",
              c(znames, ynames, xnames)
            )
          )
        )
      }
    } else {
      if (fixed_x) {
        out <- list(
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_y_alg,
            name = "expected_mean",
            dimnames = list(
              "1",
              ynames
            )
          )
        )
      } else {
        out <- list(
          OpenMx::mxAlgebraFromString(
            algString = "t(mu_x)",
            name = "expected_mean_x",
            dimnames = list(
              "1",
              xnames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = expected_mean_y_alg,
            name = "expected_mean_y",
            dimnames = list(
              "1",
              ynames
            )
          ),
          OpenMx::mxAlgebraFromString(
            algString = "cbind(expected_mean_y, expected_mean_x)",
            name = "expected_mean",
            dimnames = list(
              "1",
              c(ynames, xnames)
            )
          )
        )
      }
    }
  } else {
    expected_mean_y_alg <- "t(alpha)"
    if (distal) {
      out <- list(
        OpenMx::mxAlgebraFromString(
          algString = "t(kappa + phi %*% alpha)",
          name = "expected_mean_z",
          dimnames = list(
            "1",
            znames
          )
        ),
        OpenMx::mxAlgebraFromString(
          algString = expected_mean_y_alg,
          name = "expected_mean_y",
          dimnames = list(
            "1",
            ynames
          )
        ),
        OpenMx::mxAlgebraFromString(
          algString = "cbind(expected_mean_z, expected_mean_y)",
          name = "expected_mean",
          dimnames = list(
            "1",
            c(znames, ynames)
          )
        )
      )
    } else {
      out <- list(
        OpenMx::mxAlgebraFromString(
          algString = expected_mean_y_alg,
          name = "expected_mean",
          dimnames = list(
            "1",
            ynames
          )
        )
      )
    }
  }
  out
}
