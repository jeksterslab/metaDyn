.VHatMultivariate <- function(v,
                              p,
                              n) {
  list(
    OpenMx::mxMatrix(
      type = "Symm",
      ncol = p,
      nrow = p,
      free = FALSE,
      values = .AverageWithinMultivariate(
        v = v,
        n = n
      ),
      name = "v_hat"
    )
  )
}
