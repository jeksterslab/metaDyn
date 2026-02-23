.MxHelperSigmaDiagFromLDLSoftplusAlgebra <- function(column, # nolint: object_name_linter, line_length_linter
                                                     name) {
  OpenMx::mxAlgebraFromString(
    algString = paste0(
      "vec2diag(",
      "((",
      column,
      " + abs(",
      column,
      ")) / 2) + ",
      "mxRobustLog(1 + exp(-abs(",
      column,
      ")))",
      " + 1e-8",
      ")"
    ),
    name = name
  )
}
