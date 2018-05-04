# Set arguments for sampling
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The stanfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_*} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param prior Prior distribution list (can be NULL).
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for
#   \code{do.call(sampling, args)}.
set_sampling_args <- function(object, prior, user_dots = list(),
                              user_adapt_delta = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  defaults <- default_stan_control(prior = prior, adapt_delta = user_adapt_delta)

  if (!"control" %in% unms) {
    # no user-specified 'control' argument
    args$control <- defaults
  } else {
    # user specifies a 'control' argument
    if (!is.null(user_adapt_delta)) {
      # if user specified adapt_delta argument to stan_* then
      # set control$adapt_delta to user-specified value
      args$control$adapt_delta <- user_adapt_delta
    } else {
      # use default adapt_delta for the user's chosen prior
      args$control$adapt_delta <- defaults$adapt_delta
    }
    if (is.null(args$control$max_treedepth)) {
      # if user's 'control' has no max_treedepth set it to rstanarm default
      args$control$max_treedepth <- defaults$max_treedepth
    }
  }
  args$save_warmup <- FALSE

  return(args)
}

# Default control arguments for sampling
#
# Called by set_sampling_args to set the default 'control' argument for
# \code{rstan::sampling} if none specified by user. This allows the value of
# \code{adapt_delta} to depend on the prior.
#
# @param prior Prior distribution list (can be NULL).
# @param adapt_delta User's \code{adapt_delta} argument.
# @param max_treedepth Default for \code{max_treedepth}.
# @return A list with \code{adapt_delta} and \code{max_treedepth}.
default_stan_control <- function(prior, adapt_delta = NULL,
                                 max_treedepth = 15L) {
  if (!length(prior)) {
    if (is.null(adapt_delta)) adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist,
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "lasso" = 0.99,
                          "product_normal" = 0.99,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}

# Check whether a vector/matrix/array contains an "(Intercept)"
check_for_intercept <- function(x, logical = FALSE) {
  nms <- if (is.matrix(x)) colnames(x) else names(x)
  sel <- which("(Intercept)" %in% nms)
  if (logical) as.logical(length(sel)) else sel
}

# Drop intercept from a vector/matrix/array of named coefficients
drop_intercept <- function(x) {
  sel <- check_for_intercept(x)
  if (length(sel) && is.matrix(x)) {
    x[, -sel, drop = FALSE]
  } else if (length(sel)) {
    x[-sel]
  } else {
    x
  }
}

# Return intercept from a vector/matrix/array of named coefficients
return_intercept <- function(x) {
  sel <- which("(Intercept)" %in% names(x))
  if (length(sel)) x[sel] else NULL
}

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
  #x <- deparse(x, 500L)
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Reformulate as RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_rhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(~ RHS, list(RHS = x)))
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

# Replace an NA object, or NA entries in a vector
#
# @param x The vector with elements to potentially replace.
# @param replace_with The replacement value.
replace_na <- function(x, replace_with = "0") {
  if (is.na(x)) {
    x <- replace_with
  } else {
    x[is.na(x)] <- replace_with
  }
  x
}

# Replace an NULL object, or NULL entries in a vector
#
# @param x The vector with elements to potentially replace.
# @param replace_with The replacement value.
replace_null <- function(x, replace_with = "0") {
  if (is.null(x)) {
    x <- replace_with
  } else {
    x[is.null(x)] <- replace_with
  }
  x
}

# Replace named elements of 'x' with 'y'
replace_named_elements <- function(x, y) {
  x[names(y)] <- y
  x
}

# Check if all elements of a vector are zeros
all_zero <- function(x) {
  all(x == 0)
}

# Shorthand for as.integer, as.double, as.matrix, as.array
ai <- function(x, ...) as.integer(x, ...)
ad <- function(x, ...) as.double(x, ...)
am <- function(x, ...) as.matrix(x, ...)
aa <- function(x, ...) as.array(x, ...)

# Return a vector of 0's
zeros <- function(n) {
  rep(0, times = n)
}

# Return a vector of 1's
ones <- function(n) {
  rep(1, times = n)
}

# Return the maximum integer
max_integer <- function() {
  .Machine$integer.max
}

# Return the maximum double
max_double <- function() {
  .Machine$double.xmax
}

# Paste items, collapsed together using a comma
comma <- function(x) {
  paste(x, collapse = ", ")
}

# Check if object is a list (ensuring FALSE for data frames)
is_list <- function(x) {
  is(x, "list")
}

# Shorthand for %in%
is_in <- function(x, y) {
  x %in% y
}

# Error without printing call
stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# Warning without printing call
warning2 <- function(...) {
  warning(..., call. = FALSE)
}



