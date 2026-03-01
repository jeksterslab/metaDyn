## ---- test-metaDyn-softplus
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    if (!identical(Sys.getenv("NOT_CRAN"), "true") && !interactive()) {
      message("CRAN: tests skipped.")
      # nolint start
      return(invisible(NULL))
      # nolint end
    }
    set.seed(42)
    testthat::test_that(
      paste(text, "uc_d softplus-transforms back to d"),
      {
        testthat::skip_on_cran()
        A <- matrix(rnorm(9), 3, 3)
        X <- crossprod(A) + diag(1e-6, 3)
        out <- metaDyn:::.MxHelperLDL(X)
        testthat::expect_equal(
          metaDyn:::.MxHelperSoftplus(out$uc_d),
          out$d,
          tolerance = 1e-10
        )
      }
    )
    testthat::test_that(
      paste(text, "Softplus and InvSoftplus are numerical inverses"),
      {
        testthat::skip_on_cran()
        x <- rnorm(100)
        y <- metaDyn:::.MxHelperSoftplus(x)
        x_rec <- metaDyn:::.MxHelperInvSoftplus(y)
        testthat::expect_equal(
          x_rec,
          x,
          tolerance = 1e-10
        )
      }
    )
    testthat::test_that(
      paste(text, "Softplus output is strictly positive"),
      {
        testthat::skip_on_cran()
        x <- c(-50, -5, 0, 5, 50)
        y <- metaDyn:::.MxHelperSoftplus(x)
        testthat::expect_true(
          all(y > 0)
        )
        testthat::expect_true(
          all(is.finite(y))
        )
      }
    )
    testthat::test_that(
      paste(text, "Expect error"),
      {
        testthat::skip_on_cran()
        x <- c(-50, -5, 0, 5, 50)
        testthat::expect_error(
          metaDyn:::.MxHelperInvSoftplus(x)
        )
      }
    )
  },
  text = "test-metaDyn-softplus"
)
