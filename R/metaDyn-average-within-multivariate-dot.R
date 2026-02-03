.AverageWithinMultivariate <- function(v,
                                       n) {
  v_inv <- lapply(
    X = v,
    FUN = function(i) {
      out <- solve(i)
      0.5 * (out + t(out))
    }
  )
  sum_v_inv <- (
    1 / n
  ) * Reduce(
    f = `+`,
    x = v_inv
  )
  out <- solve(
    sum_v_inv
  )
  0.5 * (out + t(out))
}
