.MxHelperSdiagValues <- function(p,
                                 val) {
  if (is.null(val)) {
    val <- matrix(
      data = 0,
      nrow = p,
      ncol = p
    )
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
  val[
    upper.tri(
      x = val,
      diag = TRUE
    )
  ] <- 0
  val
}
