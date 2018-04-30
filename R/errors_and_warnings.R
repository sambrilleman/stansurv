STOP_named_list <- function(x) {
  if (!is_list(x)) {
    arg <- deparse(substitute(x))
    msg <- paste0("'", arg, "' should be a named list.")
    stop2(msg)
  }
}

STOP_invalid_names <- function(x, ok_nms) {
  nms <- names(x)
  if (length(nms) && !all(nms %in% ok_nms)) {
    arg <- deparse(substitute(x))
    msg <- paste0("Invalid named element specified in '", arg, "'. ",
                  "The following names are allowed: ")
    stop2(msg, commas(ok_nms))
  }
}
