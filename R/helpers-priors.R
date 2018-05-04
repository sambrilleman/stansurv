# Check the user-specified 'priors' list has valid named elements
validate_user_prior <- function(x, surv_dist = "weibull") {
  if (is.null(x)) {
    # priors was set to NULL, return empty list
    x <- list()
  } else if (!is_list(x)) {
    # priors was not a list, return error
    STOP_named_list(x)
  } else {
    # check names of the user-specified list of priors
    ok_nms <- ok_priors_parnames(dist = surv_dist)
    STOP_invalid_names(x, ok_nms)
  }
}

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
get_prior <- function(prior, nvars, default_scale,
                      ok_dists = ok_priors_for_beta()) {

  # build NULL prior if necessary
  if (!length(prior))
    return(make_null_prior(nvars))

  # check prior is a named list
  if (!is.list(prior))
    STOP_named_list(prior)

  # validate user inputs for: dist, location, scale, df
  prior_name  <- set_prior_name(prior$dist, ok_dists = ok_dists) # character
  prior_dist  <- set_prior_dist(prior$dist)                      # integer
  prior_mean  <- set_prior_mean(prior$location, nvars = nvars)
  prior_scale <- set_prior_scale(prior$scale, nvars = nvars, default_scale = default_scale)
  prior_df    <- set_prior_df(prior$df, nvars = nvars)

  # additional info for shrinkage priors
  if (prior_dist_name %in% shrink_priors()) {
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

# Return the default priors
#
# @param basehaz Character string, the distribution for the baseline hazard.
get_default_priors <- function(dist) {

  priors <- list()

  if (dist == "exponential") {
    priors$exp_scale <- exponential(rate = 1)
  } else if (dist == "weibull") {
    priors$wei_scale <- exponential(rate = 1)
    priors$wei_shape <- exponential(rate = 1)
  } else if (dist == "fpm") {
    priors$fpm_coefs <- normal(location = 0, scale = 2)
  }

  priors$betas <-

  priors
}

# Return the default scale hyperparameter for 'prior_aux'
#
# @param basehaz Character string, the distribution for the baseline hazard.
get_default_scale <- function(basehaz) {
  if (basehaz == "weibull") 2 else 20
}

# Return logical indicating whether the auxiliary parameters for the
# survival distribution has a lower bound at zero.
#
# @return TRUE for the following distributions:
#   weibull: shape parameter
bounded_aux <- function(dist) {
  dist %in% c("weibull")
}

# Return the validated prior distribution
#
# @param dist_name The user specified prior distribution.
# @param ok_dists A list of valid prior distributions.
set_prior_name(dist_name, ok_dists) {
  if (!dist_name %in% unlist(ok_dists)) {
    nms <- paste(names(ok_dists), collapse = ", ")
    stop2("The prior distribution should be one of: ", nms)
  }
  dist_name
}

# Switch prior distribution for the integer used internally by Stan
#
# @param dist Character string specifying the prior distribution.
# @return An integer, or NA if unmatched.
set_prior_dist <- function(dist_name) {
  switch(dist_name,
         normal      = 1L,
         t           = 2L,
         hs          = 3L,
         hs_plus     = 4L,
         laplace     = 5L,
         lasso       = 6L,
         exponential = 3L, # only used for scale pars, so no conflict with 3 for hs
         NA)
}

# Return the validated prior mean
#
# @param prior_mean Scalar or vector, the user specified location for the prior.
# @param nvars Integer The required length for the mean vector.
set_prior_mean <- function(prior_mean, nvars) {
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
set_prior_scale <- function(scale, nvars, default_scale) {
  stopifnot(is.numeric(default))
  scale <- replace_null(scale, replace_with = default_scale)
  scale <- maybe_broadcast(scale, nvars)
  scale
}

# Return the validated prior df
#
# @param location Scalar or vector, the user specified df for the prior.
# @param nvars Integer The required length for the df vector.
set_prior_df <- function(prior_df, nvars) {
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
ok_priors_for_beta <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy",
        "hs",
        "hs_plus",
        "laplace",
        "lasso")
}

# Return a list of valid priors for constrained auxiliary parameters
ok_priors_for_caux <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy",
        "exponential")
}

# Return a list of valid priors for unconstrained auxiliary parameters
ok_priors_for_uaux <- function(...) {
  nlist("normal",
        student_t = "t",
        "cauchy")
}

# Return the parameter names that can be specified by the user for each
# element of the 'priors' list
ok_priors_parnames <- function(dist = "exponential") {
  if (dist == "exponential") {
    out <- c("betas", "exp_scale")
  } else if (dist == "weibull") {
    out <- c("betas", "wei_shape", "wei_scale")
  } else if (dist == "fpm") {
    out <- c("betas", "fpm_coefs")
  }
  out
}

# Return the names of valid non-shrinkage priors for (unconstrained) coefficients
non_shrink_priors <- function(...) {
  setdiff(unlist(ok_priors_for_beta), shrink_priors)
}

# Return the names of valid shrinkage priors for (unconstrained) coefficients
shrink_priors <- function(...) {
  c("hs", "hs_plus")
}
