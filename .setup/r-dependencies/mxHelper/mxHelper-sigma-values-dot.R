.MxHelperSigmaValues <- function(p,
                                 val) {
  if (is.null(val)) {
    val <- diag(p)
  } else {
    if (is.vector(val)) {
      val <- matrix(
        data = val,
        nrow = p,
        ncol = p
      )
    }
    stopifnot(
      is.matrix(val),
      dim(val) == c(p, p)
    )
  }
  if (
    any(
      eigen(
        x = val,
        symmetric = TRUE,
        only.values = TRUE
      )$values <= 1e-06
    )
  ) {
    val <- as.matrix(
      Matrix::nearPD(val)$mat
    )
  }
  val
}
