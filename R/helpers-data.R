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
# @param basehaz Character string specifying the baseline hazard distribution.
# @return An integer, or NA if unmatched.
basehaz_for_stan <- function(basehaz) {
  switch(basehaz,
         exponential = 1L,
         weibull     = 2L,
         fpm         = 3L,
         NA)
}

# Deal with the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots
# @param ok_basehaz A list of admissible baseline hazards
# @param eventtime A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @return A named list with the following elements:
#   type: integer specifying the type of baseline hazard, 1L = weibull,
#     2L = b-splines, 3L = piecewise.
#   type_name: character string specifying the type of baseline hazard.
#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs).
handle_basehaz <- function(basehaz, basehaz_ops,
                           ok_basehaz = c("exponential", "weibull", "fpm"),
                           df, iknots, bknots, t_beg, t_end, d) {



  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' must be one of ", comma(ok_basehaz))

  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))

  name <- basehaz
  type <- basehaz_for_stan(basehaz)

  if (name %in% c("exponential", "weibull")) {

    if (!is.null(df))
      warning2("'df' is ignored for ", name, " baseline hazard.")
    if (!is.null(knots))
      warning2("'knots' is ignored for ", name, " baseline hazard.")

    df <- ifelse(name == "weibull", 1L, 0L) # weibull shape parameter
    user_df      <- NULL
    iknots       <- NULL
    bknots       <- NULL
    spline_type  <- NULL
    spline_basis <- NULL

  } else if (name == "fpm") {

    # uncensored event times
    tt <- t_end[d == 1]

    # knot locations for fpm baseline hazard
    degree <- 3L
    knots  <- get_knots(tt, df = df - degree, iknots = iknots, bknots = bknots)
    iknots  <- knots$iknots # internal knot locations
    bknots  <- knots$bknots # boundary knot locations

    spline_type <- "I-splines"
    spline_basis <- splines2::iSpline(tt, knots = iknots, Boundary.knots = bknots)
    df <- ncol(spline_basis) # validated df = num of basis terms
    user_df <- ncol(spline_basis)

  }

  nlist(name, type, user_df, df, iknots, bknots, spline_type, spline_basis)
}

# Return the design matrix for the baseline hazard (or log cumulative
# baseline hazard in the case of the fpm model)
#
# @param t The vector of times at which to evaluate the design matrix.
# @param basehaz A list with info about the baseline hazard, returned by a
#   call to 'handle_basehaz'
# @return A matrix
make_basehaz_x <- function(t, basehaz) {

  if (name == "weibull") {
    x <- matrix(0L, nrow = length(t), ncol = 0) # dud matrix for Stan
  } else if (name == "weibull") {
    x <- matrix(log(t), nrow = length(t), ncol = 1)
  } else if (name == "fpm") {
    basis <- basehaz$spline_basis
    if (is.null(basis))
      stop2("Bug found: could not find spline basis in 'basehaz' object.")
    x <- aa(predict(basis, t))
  }
  x
}

# Return the response vector (time)
#
# @param formula The parsed model formula.
# @param data The model frame.
# @param type The type of time variable to return.
# @return A numeric vector
make_t <- function(formula, data, type = c("beg", "end", "gap")) {

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
make_d <- function(formula, data) {

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
make_x <- function(formula, data) {

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
    if (standata$basehaz == 1) "exp_scale",
    if (standata$basehaz == 2) "wei_shape",
    if (standata$basehaz == 2) "wei_scale",
    if (standata$basehaz == 3) "fpm_coefs")
}



