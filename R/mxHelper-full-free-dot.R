#' @param m Positive integer.
#'   Number of rows.
#' @param n Positive integer.
#'   Number of columns.
#' @param val Vector or matrix with logical values.
#'   If `val = NULL`,
#'   the function will assume that
#'   all elements of the full matrix
#'   are free parameters.
#'   If `is.vector(val) == TRUE`,
#'   the function will coerce the vector
#'   into an `m` by `n` matrix.
#' @noRd
.MxHelperFullFree <- function(m,
                              n,
                              val) {
  if (is.null(val)) {
    val <- matrix(
      data = TRUE,
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
