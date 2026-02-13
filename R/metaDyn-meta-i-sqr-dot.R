.MetaISqr <- function(univariate,
                      p,
                      alpha) {
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
        name = "i_sqr_l"
      ),
      OpenMx::mxAlgebraFromString(
        algString = "solve(t(i_sqr_l)) %*% tau_sqr %*% solve(i_sqr_l)",
        name = "i_sqr_h"
      ),
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "cvectorize(diag2vec(i_sqr_h)",
          "/",
          "(diag2vec(i_sqr_h) + 1))"
        ),
        name = "i_sqr"
      )
    )
  }
  idx <- seq_len(p)
  c(
    i_sqr,
    list(
      OpenMx::mxAlgebraFromString(
        algString = "cvectorize(i_sqr)",
        name = "i_sqr_vec",
        dimnames = list(
          paste0(
            "i_sqr[",
            idx,
            ",",
            1,
            "]"
          ),
          "i_sqr_vec"
        )
      ),
      OpenMx::mxCI(
        reference = "i_sqr",
        interval = 1 - alpha
      )
    )
  )
}
