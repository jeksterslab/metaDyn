#' Summary Method for Object of Class `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param object an object of class `metadynmeta`.
#' @param alpha Numeric vector.
#'   Significance level \eqn{\alpha}.
#'   If `NULL`, the function will check `object`
#'   for `alpha` used in model fitting.
#' @param ci_type Character string.
#'   Valid values are `"wald"` and `"mc"`.
#' @param nrep Positive integer.
#'   Number of replications for `ci_type = "mc"`.
#' @param seed Random seed for `ci_type = "mc"`.
#' @param digits Integer indicating the number of decimal places to display.
#' @param robust Logical.
#'   If `TRUE`, use robust (sandwich) sampling variance-covariance matrix.
#'   If `FALSE`, use normal theory sampling variance-covariance matrix.
#'   If `NULL`, the function will check `object`
#'   if robust standard errors are available.
#' @param ... further arguments.
#'
#' @return Returns a matrix of
#'   estimates,
#'   standard errors,
#'   test statistics,
#'   degrees of freedom,
#'   p-values,
#'   and
#'   confidence intervals.
#'
#' @method summary metadynmeta
#' @keywords methods
#' @export
summary.metadynmeta <- function(object,
                                alpha = NULL,
                                ci_type = "wald",
                                robust = NULL,
                                nrep = 20000L,
                                seed = NULL,
                                digits = 4,
                                ...) {
  code <- .CheckStatusCode(
    model = object$output
  )
  stopifnot(
    ci_type %in% c("wald", "mc")
  )
  if (is.null(alpha)) {
    alpha <- object$args$alpha
  }
  if (is.null(robust)) {
    if (is.null(object$robust)) {
      robust <- FALSE
    } else {
      robust <- TRUE
    }
  }
  if (robust) {
    if (is.null(object$robust)) {
      utils::capture.output(
        suppressMessages(
          suppressWarnings(
            sandwich <- OpenMx::imxRobustSE(
              model = object$output,
              details = TRUE
            )
          )
        )
      )
    } else {
      sandwich <- object$robust
    }
    object$output@output$vcov <- sandwich$cov
    object$output@output$standardErrors <- sandwich$SE
    type <- "robust"
  } else {
    type <- "normal"
  }
  ci <- switch(ci_type,
    wald = .CIWaldMeta(
      object = object,
      alpha = alpha
    ),
    mc = .CIMCMeta(
      object = object,
      alpha = alpha,
      nrep = nrep,
      seed = seed
    )
  )
  print_summary <- round(
    x = ci,
    digits = digits
  )
  class(ci) <- c(
    "summary.metadynmeta",
    class(ci)
  )
  attr(ci, "fit") <- object
  attr(ci, "alpha") <- alpha
  attr(ci, "digits") <- digits
  attr(ci, "type") <- type
  attr(ci, "ci_type") <- ci_type
  attr(ci, "nrep") <- nrep
  attr(ci, "seed") <- seed
  attr(ci, "code") <- code
  attr(ci, "print_summary") <- print_summary
  ci
}

#' @noRd
#' @keywords internal
#' @exportS3Method print summary.metadynmeta
print.summary.metadynmeta <- function(x,
                                      ...) {
  print_summary <- attr(
    x = x,
    which = "print_summary"
  )
  object <- attr(
    x = x,
    which = "fit"
  )
  type <- attr(
    x = x,
    which = "type"
  )
  ci_type <- attr(
    x = x,
    which = "ci_type"
  )
  nrep <- attr(
    x = x,
    which = "nrep"
  )
  code <- attr(
    x = x,
    which = "code"
  )
  cat("Call:\n")
  base::print(object$call)
  cat(
    paste0(
      "\n",
      "Status code:\n",
      code,
      "\n"
    )
  )
  if (ci_type == "wald") {
    cat(
      paste0(
        "\n",
        "Wald CI ",
        "type:\n",
        "\"",
        type,
        "\"",
        "\n\n"
      )
    )
  }
  if (ci_type == "mc") {
    cat(
      paste0(
        "\n",
        "Monte Carlo CI ",
        "type:\n",
        "\"",
        type,
        "\"",
        "\n\n"
      )
    )
  }
  print(print_summary)
  invisible(object)
}

#' Print Method for Object of Class `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x an object of class `metadynmeta`.
#' @inheritParams summary.metadynmeta
#'
#' @inherit summary.metadynmeta return
#'
#' @method print metadynmeta
#' @keywords methods
#' @export
print.metadynmeta <- function(x,
                              alpha = NULL,
                              ci_type = "wald",
                              robust = NULL,
                              nrep = 20000L,
                              seed = NULL,
                              digits = 4,
                              ...) {
  stopifnot(
    ci_type %in% c("wald", "mc")
  )
  print.summary.metadynmeta(
    summary.metadynmeta(
      object = x,
      alpha = alpha,
      ci_type = ci_type,
      robust = robust,
      nrep = nrep,
      seed = seed,
      digits = digits
    )
  )
}

#' Estimated Parameter Method for an Object of Class
#' `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a vector of estimated parameters.
#'
#' @inheritParams summary.metadynmeta
#'
#' @method coef metadynmeta
#' @keywords methods
#' @export
coef.metadynmeta <- function(object,
                             ...) {
  coef(object$output)
}

#' Variance-Covariance Matrix Method for an Object of Class
#' `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns the sampling variance-covariance matrix
#'   of the estimated parameters.
#'
#' @inheritParams summary.metadynmeta
#'
#' @method vcov metadynmeta
#' @keywords methods
#' @export
vcov.metadynmeta <- function(object,
                             robust = NULL,
                             ...) {
  if (is.null(robust)) {
    if (is.null(object$robust)) {
      robust <- FALSE
    } else {
      robust <- TRUE
    }
  }
  if (robust) {
    if (is.null(object$robust)) {
      utils::capture.output(
        suppressMessages(
          suppressWarnings(
            sandwich <- OpenMx::imxRobustSE(
              model = object$output,
              details = TRUE
            )
          )
        )
      )
    } else {
      sandwich <- object$robust
    }
    object$output@output$vcov <- sandwich$cov
    object$output@output$standardErrors <- sandwich$SE
  }
  vcov(object$output)
}

#' Confidence Intervals for the Parameter Estimates
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams summary.metadynmeta
#' @param parm a specification of which parameters
#'   are to be given confidence intervals,
#'   either a vector of numbers or a vector of names.
#'   If missing, all parameters are considered.
#' @param level the confidence level required.
#' @return Returns a matrix of confidence intervals.
#'
#' @method confint metadynmeta
#' @keywords methods
#' @export
confint.metadynmeta <- function(object,
                                parm = NULL,
                                level = 0.95,
                                ci_type = "wald",
                                robust = NULL,
                                nrep = 20000L,
                                seed = NULL,
                                ...) {
  stopifnot(
    length(level) == 1
  )
  stopifnot(
    ci_type %in% c("wald", "mc")
  )
  if (is.null(robust)) {
    if (is.null(object$robust)) {
      robust <- FALSE
    } else {
      robust <- TRUE
    }
  }
  if (robust) {
    if (is.null(object$robust)) {
      utils::capture.output(
        suppressMessages(
          suppressWarnings(
            sandwich <- OpenMx::imxRobustSE(
              model = object$output,
              details = TRUE
            )
          )
        )
      )
    } else {
      sandwich <- object$robust
    }
    object$output@output$vcov <- sandwich$cov
    object$output@output$standardErrors <- sandwich$SE
  }
  ci <- switch(ci_type,
    wald = .CIWaldMeta(
      object = object,
      alpha = 1 - level
    )[, 5:6, drop = FALSE],
    mc = .CIMCMeta(
      object = object,
      alpha = 1 - level,
      nrep = nrep,
      seed = seed
    )[, 4:5, drop = FALSE]
  )
  if (is.null(parm)) {
    parameters <- rownames(
      ci
    )
    if (!is.null(parameters)) {
      parm <- parameters
    } else {
      parm <- seq_len(dim(ci)[1])
    }
  }
  ci <- ci[parm, , drop = FALSE]
  varnames <- colnames(ci)
  varnames <- gsub(
    pattern = "%",
    replacement = " %",
    x = varnames
  )
  colnames(ci) <- varnames
  ci
}

#' Extract Generic Function
#'
#' A generic function for extracting elements from objects.
#'
#' @param object An object.
#' @param what Character string.
#' @return A value determined by the specific method for the object's class.
#' @keywords methods
#' @export
extract <- function(object,
                    what) {
  UseMethod("extract")
}

#' Extract Method for an Object of Class
#' `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a list of estimates.
#'
#' @inheritParams summary.metadynmeta
#' @param what Character string.
#'   What specific matrix to extract.
#'   If `what = NULL`,
#'   extract all available matrices.
#'
#' @keywords methods
#' @export
#' @method extract metadynmeta
extract.metadynmeta <- function(object,
                                what = NULL) {
  out <- .MetaExtract(object = object)
  if (!is.null(what)) {
    stopifnot(
      what %in% c(
        "alpha",
        "beta",
        "gamma",
        "mu_x",
        "sigma_x",
        "tau_sqr",
        "tau_sqr_l",
        "tau_sqr_d",
        "i_sqr",
        "kappa",
        "phi",
        "psi",
        "psi_l",
        "psi_d",
        "omega"
      )
    )
    out <- out[[what]]
  }
  out
}
