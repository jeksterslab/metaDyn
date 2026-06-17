## ---- test-metaDyn-distal-mixed-effects-stochastic-x
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 alpha,
                 tau_sqr,
                 v_hat,
                 mu_x,
                 sigma_x,
                 gamma,
                 kappa,
                 phi,
                 omega,
                 psi) {
    message(text)
    set.seed(42)
    if (!identical(Sys.getenv("NOT_CRAN"), "true") && !interactive()) {
      message("CRAN: tests skipped.")
      return(invisible(NULL))
    }
    if (identical(Sys.getenv("GITHUB_TEST"), "true")) {
      ci <- TRUE
      n <- 1000
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
            c(
              MASS::mvrnorm(
                n = 1,
                mu = mu_x,
                Sigma = sigma_x
              )
            )
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
              mu = rep(0, length(kappa)),
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
          fixed_x = FALSE,
          robust = robust,
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
        extracted <- extract(fit)
        testthat::expect_false(
          is.null(extracted$mu_x)
        )
        testthat::expect_false(
          is.null(extracted$sigma_x)
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = c(extracted$mu_x),
                digits = 0
              ) - mu_x
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = extracted$sigma_x,
                digits = 1
              ) - sigma_x
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = mxEval(gamma, fit$output),
                digits = 1
              ) - gamma
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
  text = "test-metaDyn-distal-mixed-effects-stochastic-x",
  alpha = rep(x = 0.50, times = 2),
  tau_sqr = 0.50 * diag(2),
  v_hat = 0.10 * diag(2),
  mu_x = c(0.25, -0.25, 0.50),
  sigma_x = diag(c(0.50, 0.75, 1.00)),
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
