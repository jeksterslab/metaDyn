.CheckStatusCode <- function(model,
                             context = deparse(substitute(model))) {
  .Check <- function(x, msg) {
    if (is.null(x)) {
      stop(msg, call. = FALSE)
    } else {
      x
    }
  }
  code <- .Check(
    tryCatch(
      expr = model@output$status$code,
      error = function(e) {
        NULL
      }
    ),
    sprintf(
      "No OpenMx status code found in `%s` ",
      "(expected `model@output$status$code`).",
      context
    )
  )
  if (!is.numeric(code) || length(code) != 1L || is.na(code)) {
    stop(
      sprintf(
        "Status code in `%s` is not a single non-NA numeric value.",
        context
      ),
      call. = FALSE
    )
  }
  if (code != 0) {
    msg <- tryCatch(
      expr = model@output$status$details,
      error = function(e) {
        NULL
      }
    )
    if (is.null(msg) || !nzchar(msg)) {
      msg <- "No status details provided."
    }
    warning(
      sprintf(
        "OpenMx returned status code %s for `%s`: %s",
        code,
        context,
        msg
      ),
      call. = FALSE
    )
  }
  invisible(code)
}
