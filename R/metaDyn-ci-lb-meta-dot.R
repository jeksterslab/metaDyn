.CILBMeta <- function(object,
                      alpha) {
  out <- summary(object$output)$CI
  est <- out[, "estimate"]
  coef_names <- rownames(out)
  coef_names <- gsub(
    pattern = "^alpha_(\\d+)_(\\d+)$",
    replacement = "alpha[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^beta_(\\d+)_(\\d+)$",
    replacement = "beta[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^gamma_(\\d+)_(\\d+)$",
    replacement = "gamma[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^kappa_(\\d+)_(\\d+)$",
    replacement = "kappa[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^phi_(\\d+)_(\\d+)$",
    replacement = "phi[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^omega_(\\d+)_(\\d+)$",
    replacement = "omega[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^Model\\.",
    replacement = "",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^i_sqr_vec\\[(\\d+),(\\d+)\\]$",
    replacement = "i_sqr[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^direct\\[(\\d+),(\\d+)\\]$",
    replacement = "de[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^indirect\\[(\\d+),(\\d+)\\]$",
    replacement = "ie[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^total\\[(\\d+),(\\d+)\\]$",
    replacement = "te[\\1,\\2]",
    x = coef_names
  )
  coef_names <- gsub(
    pattern = "^ie_x(\\d+)_y(\\d+)_z(\\d+)\\[\\d+,\\d+\\]$",
    replacement = "ie_x\\1_y\\2_z\\3",
    x = coef_names
  )
  names(est) <- coef_names
  out <- do.call(
    what = "cbind",
    args = lapply(
      X = alpha,
      FUN = function(alpha) {
        current_level <- if (!is.null(object$output$intervals)) {
          object$output$intervals[[1]]$interval
        } else {
          NA_real_
        }
        if (
          isTRUE(
            all.equal(
              target = current_level,
              current = 1 - alpha
            )
          )
        ) {
          object_new <- object$output
        } else {
          object_new <- OpenMx::mxModel(
            model = object$output,
            remove = TRUE,
            object$output$intervals
          )
          object_new <- OpenMx::mxModel(
            model = object_new,
            OpenMx::mxCI(
              reference = names(object$output$intervals),
              interval = 1 - alpha
            )
          )
          object_new <- OpenMx::omxRunCI(
            model = object_new
          )
          out <- summary(object_new)$CI
          if (!all(c(out[, "note"]) == "")) {
            message(
              "Take note of the following:\n"
            )
            print(out[, "note", drop = FALSE])
          }
          ll <- out[, "lbound"]
          ul <- out[, "ubound"]
          out <- cbind(
            ll,
            ul
          )
          colnames(out) <- .ProbsofAlpha(alpha = alpha)
          out
        }
      }
    )
  )
  out <- out[, order(as.numeric(colnames(out)))]
  colnames(out) <- paste0(
    as.numeric(colnames(out)) * 100, "%"
  )
  out <- cbind(
    est = est,
    out
  )
  rownames(out) <- coef_names
  out[
    rownames(
      .CIWaldMeta(
        object = object,
        alpha = alpha
      )
    ), ,
    drop = FALSE
  ]
}
