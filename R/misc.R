# Return a list (or vector if unlist = TRUE) which
# contains the embedded elements in list x named y
fetch <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                  pad_length = NULL, unlist = FALSE) {
  ret <- lapply(x, `[[`, y)
  if (!is.null(z))
    ret <- lapply(ret, `[[`, z)
  if (!is.null(zz))
    ret <- lapply(ret, `[[`, zz)
  if (null_to_zero)
    ret <- lapply(ret, function(i) ifelse(is.null(i), 0L, i))
  if (!is.null(pad_length)) {
    padding <- rep(list(0L), pad_length - length(ret))
    ret <- c(ret, padding)
  }
  if (unlist) unlist(ret) else ret
}
# Wrapper for using fetch with unlist = TRUE
fetch_ <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                   pad_length = NULL) {
  fetch(x = x, y = y, z = z, zz = zz, null_to_zero = null_to_zero,
        pad_length = pad_length, unlist = TRUE)
}
# Wrapper for using fetch with unlist = TRUE and
# returning array. Also converts logical to integer.
fetch_array <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                        pad_length = NULL) {
  val <- fetch(x = x, y = y, z = z, zz = zz, null_to_zero = null_to_zero,
               pad_length = pad_length, unlist = TRUE)
  if (is.logical(val))
    val <- as.integer(val)
  as.array(val)
}

# Unlist the result from an lapply call
#
# @param X,FUN,... Same as lapply
uapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...))
}

# A refactored version of mapply with SIMPLIFY = FALSE
#
# @param FUN,... Same as mapply
# @param arg Passed to MoreArgs
xapply <- function(..., FUN, args = NULL) {
  mapply(FUN, ..., MoreArgs = args, SIMPLIFY = FALSE)
}

# Extract LHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
lhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[2L]]
  } else {
    out <- NULL
  }
  out
}

# Extract RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}

# Reformulate as LHS of a formula
#
# @param x A character string or expression object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_lhs <- function(x) {
  x <- deparse(x, 500L)
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Reformulate as RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_rhs <- function(x) {
  x <- deparse(x, 500L)
  x <- formula(substitute(~ RHS, list(RHS = out)))
  x
}

# Count the number of unique values
#
# @param x A vector or list
n_distinct <- function(x) {
  length(unique(x))
}

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
#
# @param ... Objects to include in the list.
# @return A named list.
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# Maybe broadcast
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make.
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

# Check if object is a list (ensuring FALSE for data frames)
is_list <- function(x) {
  is(x, "list")
}
