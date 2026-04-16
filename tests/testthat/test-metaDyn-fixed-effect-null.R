## ---- test-metaDyn-fixed-effect-null
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n,
                 alpha,
                 v_hat,
                 tol) {
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
        y <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            epsilon <- MASS::mvrnorm(
              n = 1,
              mu = c(0, 0),
              Sigma = v[[i]]
            )
            c(
              alpha + epsilon
            )
          }
        )
        fit <- Meta(
          y = y,
          v = v,
          random = FALSE,
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
        vcovs <- vcov(fit)
        testthat::expect_true(
          all(
            abs(
              round(
                x = coefs,
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
        if (requireNamespace("metaSEM")) {
          fit <- Meta(
            y = y,
            v = v,
            random = FALSE,
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
          data <- as.data.frame(
            cbind(
              y,
              v
            )
          )
          metasem <- meta(
            y = cbind(y1, y2),
            v = cbind(y1y1, y2y1, y2y2),
            data = data,
            RE.constraints = matrix(
              data = 0,
              nrow = 2,
              ncol = 2
            )
          )
          coefs_metasem <- coef(metasem)
          vcovs_metasem <- vcov(metasem)
          testthat::expect_true(
            all(
              abs(
                coefs - coefs_metasem
              ) <= 0.001
            )
          )
          testthat::expect_true(
            all(
              abs(
                c(mxEval(alpha, fit$output)) - coefs_metasem
              ) <= 0.001
            )
          )
          testthat::expect_true(
            all(
              abs(
                vcovs - vcovs_metasem
              ) <= 0.001
            )
          )
        }
      }
    )
  },
  text = "test-metaDyn-fixed-effect-null",
  alpha = rep(x = 0.50, times = 2),
  v_hat = 0.10 * diag(2)
)
