.MxHelperFullValues <- function(m,
                                n,
                                val) {
  if (is.null(val)) {
    val <- matrix(
      data = 0,
      nrow = m,
      ncol = n
    )
  } else {
    if (is.vector(val)) {
      if (length(val) == 1) {
        val <- rep(
          x = val,
          times = m * n
        )
      }
      stopifnot(
        length(val) == m * n
      )
      val <- matrix(
        data = val,
        nrow = m,
        ncol = n
      )
    }
    stopifnot(
      is.matrix(val),
      dim(val) == c(m, n)
    )
  }
  val
}
