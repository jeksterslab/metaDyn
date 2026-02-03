#' @param m Positive integer.
#'   Number of rows.
#' @param n Positive integer.
#'   Number of columns.
#' @param values Vector or matrix of starting values.
#'   If `values = NULL`,
#'   starting values will be a zero matrix.
#' @param free_val Values for `.MxHelperFullFree`.
#' @param lbound_val Lower bound values for `.MxHelperFullBound`.
#' @param ubound_val Upper bound values for `.MxHelperFullBound`.
#' @param vec Logical.
#'   If `vec = TRUE`,
#'   create a matrix of only free parameters
#'   that will be vectorized later in `.MxHelperFullMxMatrix`.
#' @param name Character string.
#'   Name of the matrix.
#' @noRd
.MxHelperFullPrepMatrices <- function(m,
                                      n,
                                      free_val,
                                      values,
                                      lbound_val,
                                      ubound_val,
                                      vec,
                                      name) {
  values <- .MxHelperFullValues(
    m = m,
    n = n,
    val = values
  )
  free <- .MxHelperFullFree(
    m = m,
    n = n,
    val = free_val
  )
  labels <- .MxHelperFullLabels(
    m = m,
    n = n,
    name = name,
    sep = "underscore"
  )
  if (vec) {
    vec <- .MxHelperFullLabels(
      m = m,
      n = n,
      name = name,
      sep = "bracket"
    )
    vec_run <- TRUE
  } else {
    vec_run <- FALSE
  }
  lbound <- .MxHelperFullBound(
    m = m,
    n = n,
    val = lbound_val
  )
  ubound <- .MxHelperFullBound(
    m = m,
    n = n,
    val = ubound_val
  )
  for (j in seq_len(n)) {
    for (i in seq_len(m)) {
      if (!free[i, j]) {
        labels[i, j] <- NA
        lbound[i, j] <- NA
        ubound[i, j] <- NA
        if (vec_run) {
          vec[i, j] <- NA
        }
      }
    }
  }
  list(
    free = free,
    values = values,
    labels = labels,
    lbound = lbound,
    ubound = ubound,
    vec = vec
  )
}
