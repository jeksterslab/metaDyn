.MxHelperSigmaDiagFromLDLMxMatrix <- function(p, # nolint: object_name_linter, line_length_linter
                                              name,
                                              column_name,
                                              d_free,
                                              d_values,
                                              d_lbound,
                                              d_ubound,
                                              d_rows,
                                              d_cols,
                                              d_equal) {
  if (is.null(d_values)) {
    d_values <- rep(
      x = log(expm1(1)),
      times = p
    )
  } else {
    if (is.matrix(d_values)) {
      if (all(dim(d_values) == c(p, p))) {
        d_values <- diag(d_values)
      }
      if (all(dim(d_values) == c(p, 1))) {
        d_values <- c(d_values)
      }
    }
    stopifnot(
      is.vector(d_values),
      length(d_values) == p
    )
  }
  out <- list()
  column <- .MxHelperFullMxMatrix(
    m = p,
    n = 1,
    free_val = d_free,
    values = d_values,
    lbound_val = d_lbound,
    ubound_val = d_ubound,
    vec = TRUE,
    row = d_rows,
    col = d_cols,
    name = column_name
  )
  if (isTRUE(d_equal)) {
    mat <- column[[1]]
    eq_label <- paste0(
      column_name,
      "_eq"
    )
    fr <- matrix(
      data = mat@free,
      nrow = p,
      ncol = 1
    )
    lab <- matrix(
      data = mat@labels,
      nrow = p,
      ncol = 1
    )
    val <- matrix(
      data = mat@values,
      nrow = p,
      ncol = 1
    )
    lab[fr] <- eq_label
    v0 <- if (length(d_values) > 0) {
      d_values[1]
    } else {
      log(expm1(1))
    }
    val[fr] <- v0
    if (!is.null(mat@lbound)) {
      lb <- matrix(
        data = mat@lbound,
        nrow = p,
        ncol = 1
      )
      lb[fr] <- if (length(d_lbound)) {
        d_lbound[1]
      } else {
        lb[fr]
      }
      mat@lbound <- lb
    }
    if (!is.null(mat@ubound)) {
      ub <- matrix(
        data = mat@ubound,
        nrow = p,
        ncol = 1
      )
      ub[fr] <- if (length(d_ubound)) {
        d_ubound[1]
      } else {
        ub[fr]
      }
      mat@ubound <- ub
    }
    mat@labels <- lab
    mat@values <- val
    column[[1]] <- mat
  }
  if (length(column) == 1) {
    names(column) <- column_name
  } else {
    names(column) <- c(
      column_name,
      paste0(
        column_name,
        "_",
        "vec"
      )
    )
  }
  out <- c(
    out,
    column
  )
  sigma <- .MxHelperSigmaDiagFromLDLSoftplusAlgebra(
    column = column_name,
    name = name
  )
  out <- c(
    out,
    stats::setNames(
      list(sigma),
      name
    )
  )
  sigma_vech <- OpenMx::mxAlgebraFromString(
    algString = paste0(
      "vech(",
      name,
      ")"
    ),
    name = paste0(
      name,
      "_",
      "vech"
    )
  )
  out <- c(
    out,
    stats::setNames(
      list(sigma_vech),
      paste0(
        name,
        "_",
        "vech"
      )
    )
  )
  free <- matrix(
    data = FALSE,
    nrow = p,
    ncol = p
  )
  diag(free) <- c(
    column[[1]]@free
  )
  vec <- matrix(
    data = NA_character_,
    nrow = p,
    ncol = p
  )
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      if (free[i, j]) {
        vec[i, j] <- paste0(
          name,
          "[",
          i,
          ",",
          j,
          "]"
        )
      }
    }
  }
  vec <- unique(
    stats::na.omit(
      c(vec)
    )
  )
  q <- length(vec)
  if (q > 0) {
    vec_free <- OpenMx::mxMatrix(
      type = "Full",
      nrow = q,
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
        list(vec_free),
        paste0(
          name,
          "_vec"
        )
      )
    )
  }
  out
}
