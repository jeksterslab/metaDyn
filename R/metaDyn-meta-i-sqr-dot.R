.MetaISqr <- function(p) {
  i_sqr <- list(
    i_sqr = OpenMx::mxAlgebraFromString(
      algString = paste0(
        "cvectorize((diag2vec(tau_sqr)",
        "/",
        "(diag2vec(tau_sqr) + v_hat)))"
      ),
      name = "i_sqr"
    )
  )
  idx <- seq_len(p)
  c(
    i_sqr,
    list(
      i_sqr_vec = OpenMx::mxAlgebraFromString(
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
      )
    )
  )
}
