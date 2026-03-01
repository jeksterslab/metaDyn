## ---- test-metaDyn-ldl
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
      paste(text, "LDL reconstructs the input matrix"),
      {
        testthat::skip_on_cran()
        a <- matrix(
          data = rnorm(16),
          nrow = 4,
          ncol = 4
        )
        x <- crossprod(a) + diag(1e-6, 4)
        ldl <- metaDyn:::.MxHelperLDL(x)
        inv_ldl <- metaDyn:::.MxHelperInvLDL(s_l = ldl$s_l, uc_d = ldl$uc_d)
        testthat::expect_true(
          max(abs(x - inv_ldl)) < 1e-8
        )
        testthat::expect_true(
          all(
            .SymofVech(.Vech(x), k = 4) == x
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "LDL diagonal entries are positive"),
      {
        testthat::skip_on_cran()
        a <- matrix(rnorm(25), 5, 5)
        x <- crossprod(a) + diag(1e-6, 5)
        out <- metaDyn:::.MxHelperLDL(x)
        testthat::expect_true(
          all(out$d > 0)
        )
      }
    )
  },
  text = "test-metaDyn-ldl"
)
