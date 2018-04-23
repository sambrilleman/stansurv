# Check input argument is a valid type, and return as a list
#
# @param arg The user input to the argument.
# @param type Character string specifying the type of argument being checked.
#   Can be either: "formula", "data", "dist".
# @param return_list Should a list be returned? Or just a single element.
# @param validate_length The required length of the returned list.
# @return A list of formulas, data frames, or distribution names.
validate_arg <- function(arg, type, return_list = FALSE,
                         validate_length = NULL, ...) {

  nm <- deparse(substitute(arg))

  ok_inputs <- switch(type,
                      formula = "formula",
                      data = "data.frame",
                      dist = "character")

  if (inherits(arg, ok_inputs)) {
    # input type is valid, so return as a list
    arg <- list(arg)
  } else if (is_list(arg)) {
    # input type is already a list, so check each element
    check <- sapply(arg, inherits, what = ok_inputs)
    if (!all(check))
      STOP_arg(nm, ok_inputs)
  } else {
    # input type is invalid
    STOP_arg(nm, ok_inputs)
  }

  if (type == "data") {
    arg <- lapply(arg, as.data.frame)
  } else if (type == "dist") {
    ok_dists <- list(...)$ok_dists
    check <- sapply(arg, function(x) x %in% ok_dists)
    if (!all(check))
      stop2("Argument 'dist' must be one of: ", paste(ok_dists, collapse = ", "))
  }

  if (!is.null(validate_length)) {
    if (length(arg) == 1L)
      arg <- rep(arg, times = validate_length)
    if (!length(arg) == validate_length)
      stop2(nm, " is a list of the incorrect length.")
  }

  if (return_list) {
    out <- arg
  } else {
    out <- arg[[1L]]
  }
  out
}



# Parse the model formula
#
# @param formula
parse_formula <- function(formula) {

  flhs <- lhs(formula)
  frhs <- rhs(formula)

  if (!inherits(flhs, "Surv")) {
    stop2("LHS of formula must be a 'Surv' object.")
  }

  nlist(lhs = flhs, rhs = frhs, all_vars = x)
}


# Return the response variable (time)
#
# @param formula The parsed model formula.
# @param data The model frame.
# @param type The type of time variable to return.
# @return A numeric vector
make_t_for_stan <- function(formula, data, type = c("beg", "end", "gap")) {

  type <- match.arg(type)

  if (formula$surv_type == "right") {
    t_beg <- rep(0, nrow(data))
    t_end <- data[[formula$tvar_end]]
  } else if (formula$surv_type == "counting") {
    t_beg <- data[[formula$tvar_beg]]
    t_end <- data[[formula$tvar_end]]
  } else {
    stop2("Cannot yet handle '", formula$surv_type, "' type Surv objects.")
  }

  if (type == "beg") {
    out <- t_beg
  } else if (type == "end") {
    out <- t_end
  } else if (type == "gap") {
    out <- t_end - t_beg
  }
  out
}

# Return the response variable (status indicator)
#
# @param formula The parsed model formula.
# @param data The model frame.
# @return A numeric vector
make_d_for_stan <- function(formula, data) {

  if (formula$surv_type == "right") {
    out <- data[[formula$dvar]]
  } else if (formula$surv_type == "counting") {
    out <- data[[formula$dvar]]
  } else {
    stop2("Cannot yet handle '", formula$surv_type, "' type Surv objects.")
  }
  out
}

# Return the fe predictor matrix
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and may have intercept.
#   xtemp: fe model matrix, centred and no intercept.
#   x_form: the formula for the fe model matrix.
#   x_bar: the column means of the model matrix.
#   has_intercept: logical for whether the submodel has an intercept
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x_for_stan <- function(formula, data) {

  x <- model.matrix(formula$fe_form, data)
  x <- drop_intercept(x)
  x_bar <- colMeans(x)
  xtemp <- sweep(xtemp, 2, x_bar, FUN = "-")

  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(xtemp, 2L, n_distinct) < 2)
  if (any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }

  nlist(x, x_bar, N = NROW(x), K = NCOL(x))
}

