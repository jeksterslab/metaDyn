## ---- test-metaDyn-mixed-effects
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 alpha,
                 tau_sqr,
                 v_hat,
                 gamma) {
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
      n <- 5000
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
        fit <- Meta(
          y = y,
          v = v,
          x = x,
          random = TRUE,
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
                x = coefs[grep("^gamma", names(coefs))],
                digits = 1
              ) - c(gamma)
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
        # benchmark with metaSEM
        ldl_tau_sqr <- metaDyn:::.MxHelperLDL(tau_sqr)
        tau_sqr_d <- ldl_tau_sqr$uc_d
        tau_sqr_l <- ldl_tau_sqr$s_l
        fit <- Meta(
          y = y,
          v = v,
          x = x,
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
          gamma_values = gamma,
          gamma_lbound = gamma - 10,
          gamma_ubound = gamma + 10,
          robust = robust,
          seed = 42
        )
        coefs <- coef(fit)
        summary_table <- summary(fit)
        y <- do.call(what = "rbind", args = y)
        colnames(y) <- c("y1", "y2")
        v <- do.call(
          what = "rbind",
          args = lapply(
            X = v,
            FUN = function(x) {
              x[
                lower.tri(
                  x = x,
                  diag = TRUE
                )
              ]
            }
          )
        )
        colnames(v) <- c("y1y1", "y2y1", "y2y2")
        x <- do.call(what = "rbind", args = x)
        colnames(x) <- c("x1", "x2", "x3")
        data <- as.data.frame(
          cbind(
            y,
            v,
            x
          )
        )
        metasem <- meta(
          y = cbind(y1, y2),
          v = cbind(y1y1, y2y1, y2y2),
          x = cbind(x1, x2, x3),
          data = data
        )
        coefs_metasem <- coef(metasem)
        vcovs_metasem <- vcov(metasem)
        summary_table_metasem <- summary(metasem)
        testthat::expect_true(
          all(
            abs(
              coefs[1:2] - coefs_metasem[1:2]
            ) <= 0.001
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(mxEval(alpha, fit$output)) - coefs_metasem[1:2]
            ) <= 0.001
          )
        )
        testthat::expect_true(
          all(
            abs(
              summary_table[1:11, 1] - coefs_metasem
            ) <= 0.001
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(
                mxEval(tau_sqr, fit$output)
              )[c(1, 2, 4)] - coefs_metasem[9:11]
            ) <= 0.001
          )
        )
        testthat::expect_true(
          all(
            abs(
              summary_table[1:11, 2] - sqrt(diag(vcovs_metasem))
            ) <= 0.001
          )
        )
      }
    )
  },
  text = "test-metaDyn-mixed-effects",
  alpha = rep(x = 0.50, times = 2),
  tau_sqr = 0.50 * diag(2),
  v_hat = 0.10 * diag(2),
  gamma = matrix(
    data = 0.50,
    nrow = 2,
    ncol = 3
  )
)
