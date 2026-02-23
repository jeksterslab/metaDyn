.MxHelperSoftplus <- function(x) {
  a <- abs(x)
  (x + a) / 2 + log1p(exp(-a))
}

.MxHelperInvSoftplus <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) {
    stop(
      ".MxHelperInvSoftplus() requires strictly positive input."
    )
  }
  out <- log(expm1(x))
  big <- x > 20
  out[big] <- x[big] + log1p(-exp(-x[big]))
  out
}
