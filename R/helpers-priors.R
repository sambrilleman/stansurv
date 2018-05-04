# Return a list with validated information about the user-specified
# prior distribution
#
# @param prior A list returned by one of the 'priors.R' functions.
# @param nvars Integer indicating the number of variables.
# @param default_scale Default value to use to scale if not specified by user.
# @param ok_dists A list of admissible distributions.
# @return A named list with the following arguments:
#   prior_dist: Integer specifying the prior distribution for Stan.
#   prior_mean: Scalar or vector of means for the prior.
#   prior_scale: Scalar or vector of scales for the prior.
#   prior_df: Scalar or vector of dfs for the prior.
#   prior_dist_name: Character string naming the prior.
#   global_*, slab_*: Additional parameters for shrinkage priors.
handle_prior <- function(prior, nvars, default_scale,
                         ok_dists = ok_dists()) {

  # build NULL prior if necessary
  if (!length(prior))
    return(make_null_prior(nvars))

  # check prior is a named list
  if (!is.list(prior))
    STOP_named_list(prior)

  # validate user inputs for: dist, location, scale, df
  prior_name  <- validate_prior_name(prior$dist, ok_dists = ok_dists) # character
  prior_dist  <- get_prior_dist(prior$dist)                           # integer
  prior_mean  <- get_prior_mean(prior$location, nvars = nvars)
  prior_scale <- get_prior_scale(prior$scale, nvars = nvars, default_scale = default_scale)
  prior_df    <- get_prior_df(prior$df, nvars = nvars)

  # additional info for shrinkage priors
  if (prior_name %in% shrink_priors()) {
    prior_scale <- NULL
    global_prior_scale <- prior$global_scale
    global_prior_df    <- prior$global_df
    slab_scale         <- prior$slab_scale
    slab_df            <- prior$slab_df
  } else {
    global_prior_scale <- 0
    global_prior_df    <- 0
    slab_scale         <- 0
    slab_df            <- 0
  }

  nlist(prior_dist,
        prior_mean,
        prior_scale,
        prior_df,
        prior_dist_name = prior_name,
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}

# Autoscaling of priors
#
# @param prior_stuff A named list returned by a call to handle_glm_prior
# @param response A vector containing the response variable, only required if
#   the priors are to be scaled by the standard deviation of the response (for
#   gaussian reponse variables only)
# @param predictors The predictor matrix, only required if the priors are to be
#   scaled by the range/sd of the predictors
# @param family A family object
# @param QR A logical specifying whether QR decomposition is used for the
#   predictor matrix
# @param min_prior_scale The minimum allowed for prior scales
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). Cannot
#   be NULL if autoscaling priors for the association parameters.
# @param ... Other arguments passed to make_assoc_terms. If autoscaling priors
#   for the association parameters then this should include 'parts' which
#   is a list containing the design matrices for the longitudinal submodel
#   evaluated at the quadrature points, as well as 'beta' and 'b' which are
#   the parameter values to use when constructing the linear predictor(s) in
#   make_assoc_terms.
# @return A named list with the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, response = NULL, predictors = NULL,
                            family = NULL, QR = FALSE, min_prior_scale = 1e-12,
                            assoc = NULL, ...) {
  ps <- prior_stuff

  if (!is.null(response) && is.gaussian(family)) {
    # use response variable for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ss <- sd(response)
      ps$prior_scale <- ss * ps$prior_scale
    }
  }

  if (!is.null(predictors) && !QR) {
    # use predictors for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ps$prior_scale <-
        pmax(min_prior_scale,
             ps$prior_scale / apply(predictors, 2L, get_scale_value))
    }
  }

  if (!is.null(assoc)) {
    # Evaluate mean and SD of each of the association terms that will go into
    # the linear predictor for the event submodel (as implicit "covariates").
    # (NB the approximate association terms are calculated using coefs
    # from the separate longitudinal submodels estimated using glmer).
    # The mean will be used for centering each association term.
    # The SD will be used for autoscaling the prior for each association parameter.
    if (is.null(family))
      stop("'family' cannot be NULL when autoscaling association parameters.")
    assoc_terms <- make_assoc_terms(family = family, assoc = assoc, ...)
    ps$a_xbar <- as.array(apply(assoc_terms, 2L, mean))
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      a_beta_scale <- apply(assoc_terms, 2L, get_scale_value)
      ps$prior_scale <- pmax(min_prior_scale, ps$prior_scale / a_beta_scale)
    }
  }

  ps$prior_scale <- as.array(pmin(.Machine$double.xmax, ps$prior_scale))
  ps
}

# Function to return the range or SD of the predictors, used for scaling the priors
# This is taken from an anonymous function in stan_glm.fit
#
# @param x A vector
get_scale_value <- function(x) {
  num.categories <- n_distinct(x)
  x.scale <- 1
  if (num.categories == 2) {
    x.scale <- diff(range(x))
  } else if (num.categories > 2) {
    x.scale <- sd(x)
  }
  return(x.scale)
}

# Return the name of the prior distribution after validating it
#
# @param dist_name The user specified prior distribution.
# @param ok_dists A vector of valid prior distributions.
validate_prior_name <- function(dist_name, ok_dists) {
  if (!dist_name %in% unlist(ok_dists)) {
    nms <- names(ok_dists)
    stop2("The prior distribution should be one of: ", comma(nms))
  }
  dist_name
}

# Switch prior distribution for the integer used internally by Stan
#
# @param dist Character string specifying the prior distribution.
# @return An integer, or NA if unmatched.
get_prior_dist <- function(dist_name) {
  switch(dist_name,
         normal      = 1L,
         t           = 2L,
         hs          = 3L,
         hs_plus     = 4L,
         laplace     = 5L,
         lasso       = 6L,
         exponential = 3L, # only used for aux pars, so no conflict with hs
         NA)
}

# Return the validated prior mean
#
# @param prior_mean Scalar or vector, the user specified location for the prior.
# @param nvars Integer The required length for the mean vector.
get_prior_mean <- function(prior_mean, nvars) {
  prior_mean <- replace_na(prior_mean, replace_with = 0)
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- aa(prior_mean) # return as array
  prior_mean
}

# Return the validated prior scale
#
# @param scale Scalar or vector, the user specified scale for the prior.
# @param nvars Integer The required length for the scale vector.
# @param default_scale Scalar, the default scale for the prior.
get_prior_scale <- function(scale, nvars, default_scale) {
  stopifnot(is.numeric(default_scale))
  scale <- replace_null(scale, replace_with = default_scale)
  scale <- maybe_broadcast(scale, nvars)
  scale
}

# Return the validated prior df
#
# @param location Scalar or vector, the user specified df for the prior.
# @param nvars Integer The required length for the df vector.
get_prior_df <- function(prior_df, nvars) {
  prior_df <- replace_na(prior_df, replace_with = 1)
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- pmin(max_double(), prior_df)
  prior_df <- aa(prior_df) # return as array
  prior_df
}

# Return dummy info for 'glm_prior()' when prior is NULL
make_null_prior <- function(nvars) {
  list(prior_dist = 0L,
       prior_dist_name = NA,
       prior_mean = aa(zeros(nvars)),
       prior_scale = aa(ones(nvars)),
       prior_df = aa(ones(nvars)),
       global_prior_scale = 0,
       global_prior_df = 0,
       slab_df = 0,
       lab_scale = 0,
       prior_autoscale = FALSE)
}

# Return a list of valid priors for (unconstrained) coefficients
ok_dists <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy",
        "hs",
        "hs_plus",
        "laplace",
        "lasso")
}

# Return a list of valid priors for auxiliary parameters
ok_dists_for_intercept <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy")
}

# Return a list of valid priors for auxiliary parameters
ok_dists_for_aux <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy",
        "exponential")
}

# Return the names of valid non-shrinkage priors for (unconstrained) coefficients
non_shrink_priors <- function(...) {
  setdiff(unlist(ok_priors_for_beta), shrink_priors)
}

# Return the names of valid shrinkage priors for (unconstrained) coefficients
shrink_priors <- function(...) {
  c("hs", "hs_plus")
}

# Return the default scale parameter for 'prior_aux'
#
# @param basehaz Character string, the distribution for the baseline hazard.
get_default_aux_scale <- function(basehaz) {
  if (basehaz == "weibull") 2 else 1
}
