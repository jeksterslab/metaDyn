.MxHelperEnsureGoodHessian <- function(model,
                                       tries_explore = 100,
                                       tries_local = 100,
                                       max_attempts = 10,
                                       grad_tol = 1e-2,
                                       hess_tol = 1e-8,
                                       eps = 1e-6,
                                       factor = 10,
                                       abs_bnd_tol = 1e-6,
                                       rel_bnd_tol = 1e-4,
                                       silent = FALSE) {
  run <- is.null(model$output) ||
    is.null(model$output$status) ||
    model$output$status$code != 0L
  if (!run) {
    good_fit <- .MxHelperIsGoodFit(
      x = model,
      tol = grad_tol
    )
    pd_hessian <- .MxHelperHasPdHessian(
      x = model,
      tol = hess_tol
    )
    bd <- .MxHelperAtBounds(
      x = model,
      abs_bnd_tol = abs_bnd_tol,
      rel_bnd_tol = rel_bnd_tol
    )$any
    run <- !(good_fit && pd_hessian && !bd)
  }
  if (run) {
    # ---- 1) Wide exploration (no Hessian) ----
    if (!silent) {
      if (interactive()) {
        cat(
          "\nStarting initial wide exploration...\n"
        )
      }
    }
    fit <- OpenMx::mxTryHardWideSearch(
      model = model,
      extraTries = tries_explore,
      checkHess = FALSE,
      silent = silent
    )
    fit <- OpenMx::mxTryHard(
      model = fit,
      extraTries = tries_local,
      silent = silent,
      jitterDistrib = "rnorm",
      scale = 0.05,
      checkHess = FALSE
    )
    if (!silent) {
      if (interactive()) {
        cat(
          "\nStarting Hessian computation...\n"
        )
      }
    }
    attempt <- 1L
    result <- NULL
    repeat {
      final <- OpenMx::mxRun(
        model = fit,
        silent = silent
      )

      if (
        .MxHelperIsGoodFit(
          x = final,
          tol = grad_tol
        ) && .MxHelperHasPdHessian(
          x = final,
          tol = hess_tol
        ) && !(
          .MxHelperAtBounds(
            x = final,
            abs_bnd_tol = abs_bnd_tol,
            rel_bnd_tol = rel_bnd_tol
          )$any
        )
      ) {
        result <- final
        break
      }

      if (attempt >= max_attempts) {
        warning(
          paste0(
            "Hessian still not positive-definite after retries;",
            " returning last fit."
          )
        )
        result <- final
        break
      }
      attempt <- attempt + 1L
      fit <- .MxHelperNudgeOffBounds(
        x = final,
        eps = eps
      )
      if (attempt == max_attempts) {
        fit <- .MxHelperRelaxBounds(
          x = fit,
          factor = factor
        )
        fit <- OpenMx::mxTryHardWideSearch(
          model = fit,
          extraTries = tries_explore,
          checkHess = FALSE,
          silent = silent
        )
      }
      fit <- OpenMx::mxTryHard(
        model = fit,
        extraTries = tries_local,
        silent = silent,
        jitterDistrib = "rnorm",
        scale = 0.05,
        checkHess = FALSE
      )
    }
    out <- result
  } else {
    out <- model
  }
  out
}
