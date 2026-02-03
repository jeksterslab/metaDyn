.CheckEstimates <- function(y,
                            v,
                            p,
                            ynames,
                            ncores) {
  len_ids <- length(y)
  # nocov start
  par <- FALSE
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (len_ids <= ncores) {
      ncores <- len_ids
    }
    if (ncores > 1) {
      par <- TRUE
    }
  }
  if (len_ids == 1) {
    par <- FALSE
  }
  foo <- function(x) {
    if (
      any(
        eigen(
          x = x,
          symmetric = TRUE,
          only.values = TRUE
        )$values <= 1e-06
      )
    ) {
      x <- as.matrix(
        Matrix::nearPD(x)$mat
      )
    }
    colnames(x) <- rownames(x) <- ynames
    x
  }
  # nocov end
  if (par) {
    # nocov start
    os_type <- Sys.info()["sysname"]
    if (os_type == "Darwin") {
      fork <- TRUE
    } else if (os_type == "Linux") {
      fork <- TRUE
    } else {
      fork <- FALSE
    }
    if (fork) {
      y <- parallel::mclapply(
        X = y,
        FUN = function(x) {
          names(x) <- ynames
          x
        },
        mc.cores = ncores
      )
      v <- parallel::mclapply(
        X = v,
        FUN = foo,
        mc.cores = ncores
      )
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit(
        parallel::stopCluster(cl = cl),
        add = TRUE
      )
      y <- parallel::parLapply(
        cl = cl,
        X = y,
        fun = function(x) {
          names(x) <- ynames
          x
        }
      )
      v <- parallel::parLapply(
        cl = cl,
        X = v,
        fun = foo
      )
    }
    # nocov end
  } else {
    y <- lapply(
      X = y,
      FUN = function(x) {
        names(x) <- ynames
        x
      }
    )
    v <- lapply(
      X = v,
      FUN = foo
    )
  }
  list(
    y = y,
    v = v
  )
}
