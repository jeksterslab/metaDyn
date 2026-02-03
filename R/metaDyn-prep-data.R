.PrepData <- function(y,
                      v,
                      x,
                      z,
                      ynames,
                      vnames,
                      xnames,
                      znames,
                      covariate,
                      distal) {
  # the dataframe is arranged as z, y, v, x
  y <- do.call(
    what = "rbind",
    args = y
  )
  colnames(y) <- ynames
  v <- do.call(
    what = "rbind",
    args = lapply(
      X = v,
      FUN = .Vech
    )
  )
  colnames(v) <- vnames
  data <- cbind(
    y,
    v
  )
  if (covariate) {
    x <- do.call(
      what = "rbind",
      args = x
    )
    colnames(x) <- xnames
    data <- cbind(
      data,
      x
    )
  }
  if (distal) {
    z <- do.call(
      what = "rbind",
      args = z
    )
    colnames(z) <- znames
    data <- cbind(
      z,
      data
    )
  }
  as.data.frame(data)
}
