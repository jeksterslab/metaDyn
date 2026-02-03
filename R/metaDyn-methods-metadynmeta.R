#' Summary Method for Object of Class `metadynmeta`
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param object an object of class `metadynmeta`.
#' @param alpha Numeric vector.
#'   Significance level \eqn{\alpha}.
#' @param digits Integer indicating the number of decimal places to display.
#' @param lb Logical.
#'   If `TRUE`, returns profile likelihood-based confidence intervals.
#'   If `FALSE`, returns Wald confidence intervals.
#' @param robust Logical.
#'   If `TRUE`, use robust (sandwich) sampling variance-covariance matrix.
#'   If `FALSE`, use normal theory sampling variance-covariance matrix.
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
                                alpha = 0.05,
                                lb = FALSE,
                                robust = FALSE,
                                digits = 4,
                                ...) {
  if (lb) {
    ci <- .CILBMeta(
      object = object,
      alpha = alpha
    )
    type <- "lb"
  } else {
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
    ci <- .CIWaldMeta(
      object = object,
      alpha = alpha
    )
  }
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
  print(
    .CheckStatusCode(
      model = object$output
    )
  )
  cat("Call:\n")
  base::print(object$call)
  cat(
    paste0(
      "\n",
      "CI ",
      "type = ",
      "\"",
      type,
      "\"",
      "\n"
    )
  )
  print(print_summary)
  invisible(x)
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
                              alpha = 0.05,
                              lb = FALSE,
                              robust = FALSE,
                              digits = 4,
                              ...) {
  print.summary.metadynmeta(
    summary.metadynmeta(
      object = x,
      alpha = alpha,
      digits = digits,
      lb = lb,
      robust = robust
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
                             robust = FALSE,
                             ...) {
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
                                lb = TRUE,
                                robust = FALSE,
                                ...) {
  stopifnot(
    length(level) == 1
  )
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
  wald <- .CIWaldMeta(
    object = object,
    alpha = 1 - level
  )[, 5:6, drop = FALSE]
  if (lb) {
    ci <- .CILBMeta(
      object = object,
      alpha = 1 - level
    )[, 2:3, drop = FALSE]
    ci <- ci[
      rownames(wald), ,
      drop = FALSE
    ]
  } else {
    ci <- wald
  }
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
        "tau_sqr",
        "i_sqr",
        "kappa",
        "phi",
        "psi",
        "omega"
      )
    )
    out <- out[[what]]
  }
  out
}
