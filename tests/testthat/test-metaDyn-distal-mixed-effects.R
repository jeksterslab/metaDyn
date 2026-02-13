## ---- test-metaDyn-distal-mixed-effects
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
      n <- 1000
      tol <- 0.10
      plus <- 0
    } else {
      ci <- FALSE
      n <- 500
      tol <- 0.50
      plus <- 2
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
        y <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            epsilon <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0),
              Sigma = v[[i]]
            )
            upsilon <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0),
              Sigma = tau_sqr
            )
            c(
              alpha + gamma %*% x[[i]] + upsilon + epsilon
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
              kappa + phi %*% y[[i]] + omega %*% x[[i]] + delta
            )
          }
        )
        tau_sqr_d_values <- c(-0.4327521, -0.4327521)
        psi_d_values <- c(-0.4327521, -0.4327521, -0.4327521, -0.4327521)
        fit <- Meta(
          y = y,
          v = v,
          x = x,
          z = z,
          random = TRUE,
          alpha_free = rep(
            x = TRUE,
            times = length(alpha)
          ),
          alpha_values = alpha,
          alpha_lbound = alpha - 10,
          alpha_ubound = alpha + 10,
          tau_sqr_diag = TRUE,
          tau_sqr_d_free = rep(
            x = TRUE,
            times = length(tau_sqr_d_values)
          ),
          tau_sqr_d_values = tau_sqr_d_values,
          tau_sqr_d_lbound = tau_sqr_d_values - 10,
          tau_sqr_d_ubound = tau_sqr_d_values + 10,
          i_sqr_univariate = FALSE,
          gamma_free = matrix(
            data = TRUE,
            nrow = nrow(gamma),
            ncol = ncol(gamma)
          ),
          gamma_values = gamma,
          gamma_lbound = gamma - 10,
          gamma_ubound = gamma + 10,
          kappa_free = rep(
            x = TRUE,
            times = length(kappa)
          ),
          kappa_values = kappa,
          kappa_lbound = kappa - 10,
          kappa_ubound = kappa + 10,
          phi_values = phi,
          phi_free = matrix(
            data = TRUE,
            nrow = nrow(phi),
            ncol = ncol(phi)
          ),
          phi_lbound = phi - 10,
          phi_ubound = phi + 10,
          omega_values = omega,
          omega_free = matrix(
            data = TRUE,
            nrow = nrow(omega),
            ncol = ncol(omega)
          ),
          omega_lbound = omega - 10,
          omega_ubound = omega + 10,
          psi_diag = TRUE,
          psi_d_free = rep(
            x = TRUE,
            times = length(psi_d_values)
          ),
          psi_d_values = psi_d_values,
          psi_d_lbound = psi_d_values - 10,
          psi_d_ubound = psi_d_values + 10,
          robust = TRUE,
          lb = TRUE,
          seed = 42
        )
        if (ci) {
          print(fit)
          vcov(fit)
          summary(fit)
          print(summary(fit))
          summary(fit, lb = TRUE)
          confint(fit)
          confint(fit, lb = TRUE)
          confint(fit, level = 0.90, lb = TRUE)
          extract(fit)
          vcov(fit, robust = TRUE)
          confint(fit, robust = TRUE)
          summary(fit, robust = TRUE)
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
                x = mxEval(v_hat, fit$output),
                digits = 1
              ) - v_hat
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
            ) <= tol + plus
          )
        )
        testthat::expect_true(
          all(
            abs(
              round(
                x = c(mxEval(kappa, fit$output)),
                digits = 0
              ) - kappa
            ) <= tol + plus
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
              ) - phi %*% gamma + omega
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-metaDyn-distal-mixed-effects",
  alpha = c(10, 10),
  tau_sqr = 0.50 * diag(2),
  v_hat = 0.10 * diag(2),
  gamma = matrix(
    data = 0.50,
    nrow = 2,
    ncol = 3
  ),
  kappa = c(10, 10, 10, 10),
  phi = matrix(
    data = 1,
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
