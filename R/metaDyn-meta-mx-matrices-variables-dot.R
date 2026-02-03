.MetaMxMatricesVariables <- function(x,
                                     v,
                                     p,
                                     m,
                                     n,
                                     xnames,
                                     vnames,
                                     random,
                                     covariate,
                                     v_hat_univariate) {
  if (covariate) {
    x_labels <- matrix(
      data = paste0(
        "data.",
        xnames
      ),
      nrow = m,
      ncol = 1
    )
    x <- OpenMx::mxMatrix(
      type = "Full",
      nrow = m,
      ncol = 1,
      free = FALSE,
      labels = x_labels,
      name = "x"
    )
  } else {
    x <- NULL
  }
  if (random) {
    v_hat <- .VHat(
      v = v,
      p = p,
      n = n,
      univariate = v_hat_univariate
    )
  } else {
    v_hat <- NULL
  }
  out <- c(
    v = .V(
      vnames = vnames,
      p = p
    ),
    v_hat = v_hat,
    x = x
  )
  out[
    !sapply(
      X = out,
      FUN = is.null
    )
  ]
}
