.VHatUnivariate <- function(v,
                            p,
                            n) {
  list(
    v_hat = OpenMx::mxMatrix(
      type = "Full",
      ncol = 1,
      nrow = p,
      free = FALSE,
      values = matrix(
        data = .AverageWithinUnivariate(
          v = v,
          n = n
        ),
        ncol = 1
      ),
      name = "v_hat"
    )
  )
}
