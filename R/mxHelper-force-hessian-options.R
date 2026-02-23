.MxHelperForceHessianOptions <- function(model) {
  model <- OpenMx::mxOption(
    model,
    "Calculate Hessian",
    "Yes"
  )
  model <- OpenMx::mxOption(
    model,
    "Hessian",
    "Yes"
  )
  model
}
