.MetaIndirect <- function(ynames,
                          xnames,
                          znames) {
  p <- length(ynames)
  m <- length(xnames)
  r <- length(znames)
  # Direct effect
  direct <- OpenMx::mxAlgebraFromString(
    algString = "omega",
    name = "direct",
    dimnames = list(
      znames,
      xnames
    )
  )
  direct_idx <- matrix(
    data = NA,
    nrow = r,
    ncol = m
  )
  for (j in seq_len(m)) {
    for (i in seq_len(r)) {
      direct_idx[i, j] <- paste0(
        "de[",
        i,
        ",",
        j,
        "]"
      )
    }
  }
  direct_idx <- c(direct_idx)
  direct_vec <- OpenMx::mxAlgebraFromString(
    algString = "cvectorize(direct)",
    name = "direct_vec",
    dimnames = list(
      direct_idx,
      "direct_vec"
    )
  )
  directs <- list(
    direct = direct,
    direct_vec = direct_vec
  )
  # Indirect effect across all mediators
  indirect <- OpenMx::mxAlgebraFromString(
    algString = "phi %*% gamma",
    name = "indirect",
    dimnames = list(znames, xnames)
  )
  indirect_idx <- matrix(
    data = NA,
    nrow = r,
    ncol = m
  )
  for (j in seq_len(m)) {
    for (i in seq_len(r)) {
      indirect_idx[i, j] <- paste0(
        "ie[",
        i,
        ",",
        j,
        "]"
      )
    }
  }
  indirect_idx <- c(indirect_idx)
  indirect_vec <- OpenMx::mxAlgebraFromString(
    algString = "cvectorize(indirect)",
    name = "indirect_vec",
    dimnames = list(
      indirect_idx,
      "indirect_vec"
    )
  )
  indirects <- list(
    indirect = indirect,
    indirect_vec = indirect_vec
  )
  # Total effect across all mediators
  total <- OpenMx::mxAlgebraFromString(
    algString = "omega + indirect",
    name = "total",
    dimnames = list(znames, xnames)
  )
  total_idx <- matrix(
    data = NA,
    nrow = r,
    ncol = m
  )
  for (j in seq_len(m)) {
    for (i in seq_len(r)) {
      total_idx[i, j] <- paste0(
        "te[",
        i,
        ",",
        j,
        "]"
      )
    }
  }
  total_idx <- c(total_idx)
  total_vec <- OpenMx::mxAlgebraFromString(
    algString = "cvectorize(total)",
    name = "total_vec",
    dimnames = list(
      total_idx,
      "total_vec"
    )
  )
  totals <- list(
    total = total,
    total_vec = total_vec
  )
  # Specific indirects: ie_<x>_<y>_<z> (1x1)
  grid_xyz <- expand.grid(
    z = seq_len(r),
    y = seq_len(p),
    x = seq_len(m),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  indirect_xyz <- lapply(
    X = seq_len(nrow(grid_xyz)),
    FUN = function(idx) {
      zi <- grid_xyz$z[idx]
      yi <- grid_xyz$y[idx]
      xi <- grid_xyz$x[idx]
      OpenMx::mxAlgebraFromString(
        algString = paste0(
          "phi[",
          zi,
          ",",
          yi,
          "] * gamma[",
          yi,
          ",",
          xi,
          "]"
        ),
        name = paste0(
          "ie_",
          xnames[xi],
          "_",
          ynames[yi],
          "_",
          znames[zi]
        ),
        dimnames = list(
          znames[zi],
          xnames[xi]
        )
      )
    }
  )
  c(
    directs,
    indirects,
    totals,
    indirect_xyz
  )
}
