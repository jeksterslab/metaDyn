#' Drop IDs From a varmxid Object
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param object Object of class `varmxid`.
#' @param drop Optional vector of unique IDs to drop.
#'
#' @return Returns `object` with selected entries in `object$converged`
#'   set to `FALSE`.
#'
#' @keywords internal
#' @noRd
.DropID <- function(object,
                    drop = NULL) {
  if (is.null(drop)) {
    return(object)
  }
  if (length(drop) < 1L) {
    return(object)
  }
  if (is.null(object$converged)) {
    stop(
      "`object$converged` is required when `drop` is not NULL.",
      call. = FALSE
    )
  }
  if (is.null(names(object$converged))) {
    stop(
      "`object$converged` must be a named logical vector when `drop` is not NULL.",
      call. = FALSE
    )
  }
  if (!is.logical(object$converged)) {
    stop(
      "`object$converged` must be logical.",
      call. = FALSE
    )
  }

  drop <- unique(as.character(drop))
  id_names <- names(object$converged)
  id_unique <- sub(
    pattern = "^(CTVAR_ID|DTVAR_ID)",
    replacement = "",
    x = id_names
  )

  drop_id <- id_names %in% drop | id_unique %in% drop

  if (!any(drop_id)) {
    warning(
      "No entries in `drop` matched names in `object$converged`.",
      call. = FALSE
    )
    return(object)
  }

  object$converged[drop_id] <- FALSE
  object
}
