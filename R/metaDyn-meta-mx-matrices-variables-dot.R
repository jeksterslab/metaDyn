.MetaMxMatricesVariables <- function(x,
                                     v,
                                     p,
                                     m,
                                     n,
                                     xnames,
                                     vnames,
                                     random,
                                     covariate,
                                     fixed_x) {
  if (covariate) {
    if (fixed_x) {
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
      # no need to define data.x if fixed_x = FALSE
      x <- NULL
    }
  } else {
    x <- NULL
  }
  if (random) {
    v_hat <- .VHat(
      v = v,
      p = p,
      n = n
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
