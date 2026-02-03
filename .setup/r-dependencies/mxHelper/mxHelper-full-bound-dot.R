#' @param m Positive integer.
#'   Number of rows.
#' @param n Positive integer.
#'   Number of columns.
#' @param val Vector or matrix with numeric values.
#'   If `(is.null(val) || all(is.na(val)))`
#'   the function will assume that
#'   the values are `NA`.
#'   If `is.vector(val) == TRUE`,
#'   the function will coerce the vector
#'   into an `m` by `n` matrix.
#' @noRd
.MxHelperFullBound <- function(m,
                               n,
                               val) {
  if (is.null(val) || all(is.na(val))) {
    val <- matrix(
      data = NA,
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
