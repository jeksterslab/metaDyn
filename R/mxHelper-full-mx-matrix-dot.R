#' @param m Positive integer.
#'   Number of rows.
#' @param n Positive integer.
#'   Number of columns.
#' @param values Numeric matrix.
#'   Matrix of starting values.
#' @param free_val Values for `.MxHelperFullFree`.
#' @param lbound_val Lower bound values for `.MxHelperFullBound`.
#' @param ubound_val Upper bound values for `.MxHelperFullBound`.
#' @param vec Logical.
#'   If `vec = TRUE`,
#'   create a matrix of only free parameters
#'   that will be vectorized later in `.MxHelperFullMxMatrix`.
#' @param row Character vector of row names.
#'   If `row = NULL`,
#'   the function will use `i_1`, `i_2`, \cdots `i_m`.
#' @param col Character vector of column names.
#'   If `col = NULL`,
#'   the function will use `j_1`, `j_2`, \cdots `i_n`.
#' @param name Character string.
#'   Name of the matrix.
#' @noRd
.MxHelperFullMxMatrix <- function(m,
                                  n,
                                  free_val,
                                  values,
                                  lbound_val,
                                  ubound_val,
                                  vec,
                                  row,
                                  col,
                                  name) {
  matrices <- .MxHelperFullPrepMatrices(
    m = m,
    n = n,
    free_val = free_val,
    values = values,
    lbound_val = lbound_val,
    ubound_val = ubound_val,
    vec = vec,
    name = name
  )
  if (is.null(row)) {
    row <- paste0(
      "i_",
      seq_len(m)
    )
  }
  if (is.null(col)) {
    col <- paste0(
      "j_",
      seq_len(n)
    )
  }
  out <- list()
  mat <- OpenMx::mxMatrix(
    type = "Full",
    nrow = m,
    ncol = n,
    free = matrices$free,
    values = matrices$values,
    labels = matrices$labels,
    lbound = matrices$lbound,
    ubound = matrices$ubound,
    byrow = FALSE,
    dimnames = list(
      row,
      col
    ),
    name = name
  )
  out <- c(
    out,
    stats::setNames(
      object = list(mat),
      nm = name
    )
  )
  if (is.matrix(matrices$vec)) {
    vec <- unique(
      c(
        stats::na.omit(
          c(
            matrices$vec
          )
        )
      )
    )
    p <- length(vec)
    if (p > 0) {
      vec_free <- OpenMx::mxMatrix(
        type = "Full",
        nrow = p,
        ncol = 1,
        labels = vec,
        dimnames = list(
          vec,
          paste0(
            name,
            "_vec"
          )
        ),
        name = paste0(
          name,
          "_vec"
        )
      )
      out <- c(
        out,
        stats::setNames(
          object = list(vec_free),
          nm = paste0(
            name,
            "_vec"
          )
        )
      )
    }
  }
  out
}
