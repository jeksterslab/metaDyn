.MxHelperSdiagPrepMatrices <- function(p,
                                       free_val,
                                       values,
                                       lbound_val,
                                       ubound_val,
                                       vec,
                                       name) {
  values <- .MxHelperSdiagValues(
    p = p,
    val = values
  )
  free <- .MxHelperSdiagFree(
    p = p,
    val = free_val
  )
  labels <- .MxHelperSdiagLabels(
    p = p,
    name = name,
    sep = "underscore"
  )
  if (vec) {
    vec <- .MxHelperSdiagLabels(
      p = p,
      name = name,
      sep = "bracket"
    )
    vec_run <- TRUE
  } else {
    vec_run <- FALSE
  }
  lbound <- .MxHelperSdiagBound(
    p = p,
    val = lbound_val
  )
  ubound <- .MxHelperSdiagBound(
    p = p,
    val = ubound_val
  )
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
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
