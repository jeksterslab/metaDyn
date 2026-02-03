#' @param m Positive integer.
#'   Number of rows.
#' @param n Positive integer.
#'   Number of columns.
#' @param name Character string.
#'   Name of the matrix.
#' @param sep Character string.
#'   Valid values are `"underscore"` and `"bracket"`.
#'   Separator used in naming free parameters.
#' @noRd
.MxHelperFullLabels <- function(m,
                                n,
                                name,
                                sep) {
  m_idx <- seq_len(m)
  n_idx <- seq_len(n)
  if (sep == "bracket") {
    out <- outer(
      X = m_idx,
      Y = n_idx,
      FUN = function(x, y) {
        paste0(
          name,
          "[",
          x,
          ",",
          y,
          "]"
        )
      }
    )
  }
  if (sep == "underscore") {
    out <- outer(
      X = m_idx,
      Y = n_idx,
      FUN = function(x, y) {
        paste0(
          name,
          "_",
          x,
          "_",
          y
        )
      }
    )
  }
  out
}
