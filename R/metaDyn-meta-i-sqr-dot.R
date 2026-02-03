.MetaISqr <- function(univariate) {
  if (univariate) {
    i_sqr <- list(
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "cvectorize((diag2vec(tau_sqr)",
          "/",
          "(diag2vec(tau_sqr) + v_hat)))"
        ),
        name = "i_sqr"
      )
    )
  } else {
    i_sqr <- list(
      OpenMx::mxAlgebraFromString(
        algString = "chol(v_hat)",
        name = "l"
      ),
      OpenMx::mxAlgebraFromString(
        algString = "solve(t(l)) %*% tau_sqr %*% solve(l)",
        name = "h"
      ),
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "cvectorize(diag2vec(h)",
          "/",
          "(diag2vec(h) + 1))"
        ),
        name = "i_sqr"
      )
    )
  }
  c(
    i_sqr,
    list(
      OpenMx::mxAlgebraFromString(
        algString = "cvectorize(i_sqr)",
        name = "i_sqr_vec"
      ),
      OpenMx::mxCI(
        reference = "i_sqr",
        interval = 0.95
      )
    )
  )
}
