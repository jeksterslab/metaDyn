.MxHelperLDL <- function(x,
                         epsilon = 1e-10) {
  diag(x) <- ifelse(
    test = diag(x) == 0,
    yes = epsilon,
    no = diag(x)
  )
  k <- dim(x)[1]
  l <- diag(k)
  d <- numeric(k)
  for (i in seq_len(k)) {
    # compute d_mat[i]
    if (i == 1) {
      d[i] <- x[i, i]
    } else {
      s <- 0.0
      for (p in seq_len(i - 1)) {
        s <- s + l[i, p] * l[i, p] * d[p]
      }
      d[i] <- x[i, i] - s
    }
    if (i < k) {
      for (j in (i + 1):k) {
        s <- 0.0
        if (i > 1) {
          for (p in seq_len(i - 1)) {
            s <- s + l[j, p] * l[i, p] * d[p]
          }
        }
        l[j, i] <- (x[j, i] - s) / d[i]
      }
    }
  }
  list(
    l = l,
    s_l = l - diag(k),
    d = d,
    uc_d = .MxHelperInvSoftplus(d),
    x = x,
    epsilon = epsilon
  )
}

.MxHelperInvLDL <- function(s_l,
                            uc_d) {
  iden <- d <- diag(dim(s_l)[1])
  s_l[
    upper.tri(
      x = s_l,
      diag = TRUE
    )
  ] <- 0
  l <- s_l + iden
  diag(d) <- .MxHelperSoftplus(c(uc_d))
  sigma <- l %*% d %*% t(l)
  0.5 * (sigma + t(sigma))
}
