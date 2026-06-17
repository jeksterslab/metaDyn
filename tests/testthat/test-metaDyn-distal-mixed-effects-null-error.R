## ---- test-metaDyn-distal-mixed-effects-null-error
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
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            alpha_free = rep(
              x = FALSE,
              times = 2
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            tau_sqr_d_free = rep(
              x = FALSE,
              times = 2
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            tau_sqr_l_free = rep(
              x = FALSE,
              times = 2 * 2
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            gamma_free = rep(
              x = FALSE,
              times = 2 * 3
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            phi_free = rep(
              x = FALSE,
              times = 4 * 2
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            omega_free = rep(
              x = FALSE,
              times = 4 * 3
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            psi_d_free = rep(
              x = FALSE,
              times = 4
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42,
            psi_diag = FALSE,
            psi_l_free = rep(
              x = FALSE,
              times = 4 * 4
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = FALSE,
            seed = 42,
            mu_x_free = rep(
              x = FALSE,
              times = 3
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = FALSE,
            seed = 42,
            sigma_x_d_free = rep(
              x = FALSE,
              times = 3
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = x,
            z = z,
            random = TRUE,
            fixed_x = FALSE,
            seed = 42,
            sigma_x_l_free = matrix(
              data = FALSE,
              nrow = 3,
              ncol = 3
            )
          )
        )
        testthat::expect_error(
          Meta(
            y = y,
            v = v,
            x = lapply(
              X = x,
              FUN = function(xi) {
                xi[1] <- NA_real_
                xi
              }
            ),
            z = z,
            random = TRUE,
            fixed_x = TRUE,
            seed = 42
          ),
          regexp = "missing values in x"
        )
      }
    )
  },
  text = "test-metaDyn-distal-mixed-effects-null-error",
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
