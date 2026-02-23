.MxHelperSdiagMxMatrix <- function(p,
                                   free_val,
                                   values,
                                   lbound_val,
                                   ubound_val,
                                   vec,
                                   row,
                                   col,
                                   name) {
  matrices <- .MxHelperSdiagPrepMatrices(
    p = p,
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
      seq_len(p)
    )
  }
  if (is.null(col)) {
    col <- paste0(
      "j_",
      seq_len(p)
    )
  }
  out <- list()
  mat <- OpenMx::mxMatrix(
    type = "Sdiag",
    nrow = p,
    ncol = p,
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
