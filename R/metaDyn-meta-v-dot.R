.V <- function(vnames,
               p) {
  list(
    OpenMx::mxMatrix(
      type = "Symm",
      nrow = p,
      ncol = p,
      labels = paste0(
        "data.",
        vnames
      ),
      name = "v"
    )
  )
}
