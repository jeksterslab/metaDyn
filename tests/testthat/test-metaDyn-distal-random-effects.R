## ---- test-metaDyn-distal-random-effects
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 alpha,
                 tau_sqr,
                 v_hat,
                 kappa,
                 phi,
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
      n <- 500
      tol <- 0.10
      plus <- 0
    } else {
      ci <- FALSE
      n <- 100
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
              alpha + upsilon + epsilon
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
              kappa + phi %*% y[[i]] + delta
            )
          }
        )
        tau_sqr_d_values <- c(-0.4327521, -0.4327521)
        tau_sqr_l_values <- matrix(
          data = 0,
          nrow = nrow(tau_sqr),
          ncol = ncol(tau_sqr)
        )
        psi_d_values <- c(-0.4327521, -0.4327521, -0.4327521, -0.4327521)
        psi_l_values <- matrix(
          data = 0,
          nrow = nrow(psi),
          ncol = ncol(psi)
        )
        fit <- Meta(
          y = y,
          v = v,
          z = z,
          random = TRUE,
          alpha_free = rep(
            x = TRUE,
            times = length(alpha)
          ),
          alpha_values = alpha,
          alpha_lbound = alpha - 10,
          alpha_ubound = alpha + 10,
          tau_sqr_diag = FALSE,
          tau_sqr_d_free = rep(
            x = TRUE,
            times = length(tau_sqr_d_values)
          ),
          tau_sqr_d_values = tau_sqr_d_values,
          tau_sqr_d_lbound = tau_sqr_d_values - 10,
          tau_sqr_d_ubound = tau_sqr_d_values + 10,
          tau_sqr_l_free = matrix(
            data = TRUE,
            nrow = nrow(tau_sqr),
            ncol = ncol(tau_sqr)
          ),
          tau_sqr_l_values = tau_sqr_l_values,
          tau_sqr_l_lbound = tau_sqr_l_values - 10,
          tau_sqr_l_ubound = tau_sqr_l_values + 10,
          i_sqr_univariate = FALSE,
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
          psi_diag = FALSE,
          psi_d_free = rep(
            x = TRUE,
            times = length(psi_d_values)
          ),
          psi_d_values = psi_d_values,
          psi_d_lbound = psi_d_values - 10,
          psi_d_ubound = psi_d_values + 10,
          psi_l_free = matrix(
            data = TRUE,
            nrow = nrow(psi),
            ncol = ncol(psi)
          ),
          psi_l_values = psi_l_values,
          psi_l_lbound = psi_l_values - 10,
          psi_l_ubound = psi_l_values + 10,
          robust = TRUE,
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
                x = mxEval(psi, fit$output),
                digits = 1
              ) - psi
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-metaDyn-distal-random-effects",
  alpha = c(10, 10),
  tau_sqr = 0.50 * diag(2),
  v_hat = 0.10 * diag(2),
  kappa = c(10, 10, 10, 10),
  phi = matrix(
    data = 1,
    nrow = 4,
    ncol = 2
  ),
  psi = 0.50 * diag(4)
)
