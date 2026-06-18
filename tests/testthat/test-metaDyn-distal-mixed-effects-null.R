## ---- test-metaDyn-distal-mixed-effects-null
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 alpha,
                 tau_sqr,
                 v_hat,
                 gamma,
                 kappa,
                 phi,
                 omega,
                 psi) {
    message(text)
    set.seed(42)
    if (!identical(Sys.getenv("NOT_CRAN"), "true") && !interactive()) {
      message("CRAN: tests skipped.")
      # nolint start
      return(invisible(NULL))
      # nolint end
    }
    if (identical(Sys.getenv("GITHUB_TEST"), "true")) {
      ci <- TRUE
      n <- 3000
      robust <- TRUE
      tol <- 0.50
    } else {
      ci <- FALSE
      n <- 500
      robust <- FALSE
      tol <- 0.50
    }
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        v <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            cov(
              MASS::mvrnorm(
                n = 100,
                mu = c(0, 0),
                Sigma = v_hat
              )
            )
          }
        )
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            rnorm(n = 3)
          }
        )
        eta <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            upsilon <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0),
              Sigma = tau_sqr
            )
            c(
              alpha + gamma %*% x[[i]] + upsilon
            )
          }
        )
        y <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            epsilon <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0),
              Sigma = v[[i]]
            )
            c(
              eta[[i]] + epsilon
            )
          }
        )
        z <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            delta <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0, 0, 0),
              Sigma = psi
            )
            c(
              kappa + phi %*% eta[[i]] + omega %*% x[[i]] + delta
            )
          }
        )
        fit <- Meta(
          y = y,
          v = v,
          x = x,
          z = z,
          random = TRUE,
          fixed_x = TRUE,
          seed = 42
        )
        if (ci) {
          print(fit)
          vcov(fit)
          summary(fit)
          print(summary(fit))
          confint(fit)
          extract(fit)
          vcov(fit, robust = TRUE)
          confint(fit, robust = TRUE, ci_type = "mc")
          summary(fit, robust = TRUE, ci_type = "mc")
        }
        coefs <- coef(fit)
        testthat::expect_true(
          all(
            abs(
              round(
                x = coefs[grep("^alpha", names(coefs))],
                digits = 0
              ) - alpha
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = c(mxEval(alpha, fit$output)),
                digits = 0
              ) - alpha
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(tau_sqr, fit$output),
                digits = 1
              ) - tau_sqr
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = c(
                  mxEval(v_hat, fit$output)
                ),
                digits = 1
              ) - c(
                diag(v_hat)
              )
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = coefs[grep("^kappa", names(coefs))],
                digits = 0
              ) - kappa
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = c(mxEval(kappa, fit$output)),
                digits = 0
              ) - kappa
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = coefs[grep("^phi", names(coefs))],
                digits = 0
              ) - c(phi)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(phi, fit$output),
                digits = 0
              ) - phi
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = coefs[grep("^omega", names(coefs))],
                digits = 1
              ) - c(omega)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(omega, fit$output),
                digits = 1
              ) - omega
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(psi, fit$output),
                digits = 1
              ) - psi
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(direct, fit$output),
                digits = 2
              ) - omega
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(indirect, fit$output),
                digits = 2
              ) - phi %*% gamma
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(total, fit$output),
                digits = 2
              ) - (phi %*% gamma + omega)
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-metaDyn-distal-mixed-effects-null",
  alpha = rep(x = 0.50, times = 2),
  tau_sqr = 0.50 * diag(2),
  v_hat = 0.10 * diag(2),
  gamma = matrix(
    data = 0.50,
    nrow = 2,
    ncol = 3
  ),
  kappa = rep(x = 0.50, times = 4),
  phi = matrix(
    data = 0.50,
    nrow = 4,
    ncol = 2
  ),
  omega = matrix(
    data = 0.20,
    nrow = 4,
    ncol = 3
  ),
  psi = 0.50 * diag(4)
)
