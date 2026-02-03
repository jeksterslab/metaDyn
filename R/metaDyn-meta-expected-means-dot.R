.MetaExpectedMeans <- function(covariate,
                               distal,
                               ynames,
                               znames) {
  if (covariate) {
    y1 <- "t(alpha + gamma %*% x)"
    if (distal) {
      out <- list(
        OpenMx::mxAlgebraFromString(
          algString = paste0(
            "t(",
            "kappa",
            " + ",
            "phi %*% (alpha + gamma %*% x)",
            " + ",
            "(omega %*% x)",
            ")"
          ),
          name = "expected_mean_z",
          dimnames = list(
            "1",
            znames
          )
        ),
        OpenMx::mxAlgebraFromString(
          algString = y1,
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
          algString = y1,
          name = "expected_mean",
          dimnames = list(
            "1",
            ynames
          )
        )
      )
    }
  } else {
    y1 <- "t(alpha)"
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
          algString = y1,
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
          algString = y1,
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
