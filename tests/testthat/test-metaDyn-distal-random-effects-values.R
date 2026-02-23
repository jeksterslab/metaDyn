## ---- test-metaDyn-distal-random-effects-values
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
      n <- 1000
      robust <- TRUE
      tol <- 0.10
      plus <- 0
    } else {
      ci <- FALSE
      n <- 500
      robust <- FALSE
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
        ldl_tau_sqr <- metaDyn:::.MxHelperLDL(tau_sqr)
        tau_sqr_d <- ldl_tau_sqr$uc_d
        tau_sqr_l <- ldl_tau_sqr$s_l
        ldl_psi <- metaDyn:::.MxHelperLDL(psi)
        psi_d <- ldl_psi$uc_d
        psi_l <- ldl_psi$s_l
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
          tau_sqr_d_free = rep(
            x = TRUE,
            times = length(tau_sqr_d)
          ),
          tau_sqr_d_values = tau_sqr_d,
          tau_sqr_d_lbound = -30,
          tau_sqr_d_ubound = 600,
          tau_sqr_l_free = matrix(
            data = TRUE,
            nrow = nrow(tau_sqr_l),
            ncol = ncol(tau_sqr_l)
          ),
          tau_sqr_l_values = tau_sqr_l,
          tau_sqr_l_lbound = -10,
          tau_sqr_l_ubound = 10,
          kappa_values = kappa,
          kappa_lbound = kappa - 10,
          kappa_ubound = kappa + 10,
          phi_values = phi,
          phi_lbound = phi - 10,
          phi_ubound = phi + 10,
          psi_d_values = psi_d,
          psi_d_lbound = -30,
          psi_d_ubound = 600,
          psi_l_free = matrix(
            data = TRUE,
            nrow = nrow(psi_l),
            ncol = ncol(psi_l)
          ),
          psi_l_values = psi_l,
          psi_l_lbound = -10,
          psi_l_ubound = 10,
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
  text = "test-metaDyn-distal-random-effects-values",
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
