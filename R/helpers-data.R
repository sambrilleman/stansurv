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
      stop2("Argument '", nm, "' must be one of: ", paste(ok_dists, collapse = ", "))
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
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {

  formula <- validate_formula(formula, needs_response = TRUE)

  lhs <- lhs(formula) # full LHS of formula
  rhs <- rhs(formula) # full RHS of formula

  lhs_form <- reformulate_lhs(lhs)
  rhs_form <- reformulate_rhs(rhs)

  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)

  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")

  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
  }

  nlist(lhs = lhs,
        rhs = rhs,
        lhs_form = lhs_form,
        rhs_form = rhs_form,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars = allvars,
        allvars_form = allvars_form,
        tvar_beg = tvar_beg,
        tvar_end = tvar_end,
        dvar = dvar,
        surv_type = attr(surv, "type"))
}

# Check formula object
#
# @param formula The user input to the formula argument.
# @param needs_response A logical; if TRUE then formula must contain a LHS.
validate_formula <- function(formula, needs_response = TRUE) {

  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }

  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object; the LHS of a formula evaluated in a data frame environment.
# @param ok_types A character vector giving the allowed types of Surv object.
validate_surv <- function(x, ok_types = c("right", "counting")) {

  if (!inherits(x, "Surv")) {
    stop2("LHS of 'formula' must be a 'Surv' object.")
  }

  if (!attr(x, "type") %in% ok_types) {
    stop2("Surv object type must be one of: ", comma(ok_types))
  }
  x
}

# Switch survival distribution for integer used internally by Stan
#
# @param dist Character string specifying the survival distribution.
# @return An integer, or NA if unmatched.
surv_dist_for_stan <- function(dist) {
  switch(dist,
         exponential = 1L,
         weibull     = 2L,
         fpm         = 3L,
         NA)
}

# Return the response vector (time)
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

# Return the response vector (status indicator)
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
    stop2("Bug found: cannot handle '", formula$surv_type, "' Surv objects.")
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
  xbar <- colMeans(x)

  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(x, 2L, n_distinct) < 2)
  if (any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }

  nlist(x, xbar, N = NROW(x), K = NCOL(x))
}

# Return the list of pars for Stan to monitor
#
# @param standata The list of data to pass to Stan.
# @return A character vector
pars_to_monitor <- function(standata) {
  c(if (standata$K > 0) "beta",
    if (standata$dist == 1) "exp_scale",
    if (standata$dist == 2) "wei_shape",
    if (standata$dist == 2) "wei_scale",
    if (standata$dist == 3) "fpm_coefs")
}

