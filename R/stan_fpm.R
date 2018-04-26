#' Bayesian proportional hazards regression
#'
#' Bayesian inference for proportional hazards regression models. The user
#' can specify a variety of standard parametric distributions for the
#' baseline hazard, or a Royston-Parmar flexible parametric model.
#'
#' @export
#'
#' @param formula A two-sided linear formula object describing both the
#'   fixed-effects and random-effects parts of the longitudinal submodel
#'   similar in vein to formula specification in the \strong{lme4} package
#'   (see \code{\link[lme4]{glmer}} or the \strong{lme4} vignette for details).
#'   Note however that the double bar (\code{||}) notation is not allowed
#'   when specifying the random-effects parts of the formula, and neither
#'   are nested grouping factors (e.g. \code{(1 | g1/g2))} or
#'   \code{(1 | g1:g2)}, where \code{g1}, \code{g2} are grouping factors.
#'   For a multivariate GLM this should be a list of such formula objects,
#'   with each element of the list providing the formula for one of the
#'   GLM submodels.
#' @param data A data frame containing the variables specified in
#'   \code{formula}. For a multivariate GLM, this can
#'   be either a single data frame which contains the data for all
#'   GLM submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the GLM submodels.
#' @param family The family (and possibly also the link function) for the
#'   GLM submodel(s). See \code{\link[lme4]{glmer}} for details.
#'   If fitting a multivariate GLM, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the GLM submodels. In other words, a different family
#'   can be specified for each GLM submodel.
#' @param weights Same as in \code{\link[stats]{glm}},
#'   except that when fitting a multivariate GLM and a list of data frames
#'   is provided in \code{data} then a corresponding list of weights
#'   must be provided. If weights are
#'   provided for one of the GLM submodels, then they must be provided for
#'   all GLM submodels.
#' @param prior,prior_intercept,prior_aux Same as in \code{\link{stan_glmer}}
#'   except that for a multivariate GLM a list of priors can be provided for
#'   any of \code{prior}, \code{prior_intercept} or \code{prior_aux} arguments.
#'   That is, different priors can optionally be specified for each of the GLM
#'   submodels. If a list is not provided, then the same prior distributions are
#'   used for each GLM submodel. Note that the \code{"product_normal"} prior is
#'   not allowed for \code{stan_mvmer}.
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{priors}} for
#'   more information about the prior distributions on covariance matrices.
#'   Note however that the default prior for covariance matrices in
#'   \code{stan_mvmer} is slightly different to that in \code{\link{stan_glmer}}
#'   (the details of which are described on the \code{\link{priors}} page).
#' @param init The method for generating initial values. See
#'   \code{\link[rstan]{stan}}.
#'
#' @details The \code{stan_mvmer} function can be used to fit a multivariate
#'   generalized linear model (GLM) with group-specific terms. The model consists
#'   of distinct GLM submodels, each which contains group-specific terms; within
#'   a grouping factor (for example, patient ID) the grouping-specific terms are
#'   assumed to be correlated across the different GLM submodels. It is
#'   possible to specify a different outcome type (for example a different
#'   family and/or link function) for each of the GLM submodels. \cr
#'   \cr
#'   Bayesian estimation of the model is performed via MCMC, in the same way as
#'   for \code{\link{stan_glmer}}. Also, similar to \code{\link{stan_glmer}},
#'   an unstructured covariance matrix is used for the group-specific terms
#'   within a given grouping factor, with priors on the terms of a decomposition
#'   of the covariance matrix.See \code{\link{priors}} for more information about
#'   the priors distributions that are available for the covariance matrices,
#'   the regression coefficients and the intercept and auxiliary parameters.
#'
#' @return A \link[=stanreg-objects]{stanmvreg} object is returned.
#'
#' @seealso \code{\link{stan_glmer}}, \code{\link{stan_jm}},
#'   \code{\link{stanreg-objects}}, \code{\link{stanmvreg-methods}},
#'   \code{\link{print.stanmvreg}}, \code{\link{summary.stanmvreg}},
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}}.
#'
#' @examples
#' \donttest{
#' #####
#' # A multivariate GLM with two submodels. For the grouping factor 'id', the
#' # group-specific intercept from the first submodel (logBili) is assumed to
#' # be correlated with the group-specific intercept and linear slope in the
#' # second submodel (albumin)
#' f1 <- stan_mvmer(
#'         formula = list(
#'           logBili ~ year + (1 | id),
#'           albumin ~ sex + year + (year | id)),
#'         data = pbcLong,
#'         # this next line is only to keep the example small in size!
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f1)
#'
#' #####
#' # A multivariate GLM with one bernoulli outcome and one
#' # gaussian outcome. We will artificially create the bernoulli
#' # outcome by dichotomising log serum bilirubin
#' pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
#' f2 <- stan_mvmer(
#'         formula = list(
#'           ybern ~ year + (1 | id),
#'           albumin ~ sex + year + (year | id)),
#'         data = pbcLong,
#'         family = list(binomial, gaussian),
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' }
#'
stan_fpm <- function(formula, data, dist = "fpm", df = 3, iknots = NULL, bknots = NULL,
                     priors = list(), prior_PD = FALSE,
                     algorithm = c("sampling", "meanfield", "fullrank"),
                     adapt_delta = 0.95, max_treedepth = 11L,
                     init = "random", ...) {

  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------

  dots <- list(...)
  algorithm <- match.arg(algorithm)

  # Formula
  formula <- parse_formula(formula, data)

  # Data
  data <- validate_arg(data, "data")

  # Distribution
  ok_dists <- c("exponential", "weibull", "fpm")
  dist <- validate_arg(dist, "dist", ok_dists = ok_dists)

  #----------------
  # Construct data
  #----------------

  standata <- list()

  # model data frame
  mf <- data

  # survival time distribution
  standata$dist <- make_dist_for_stan(dist)

  # time variable for each row of data
  standata$t_beg <- make_t_for_stan(formula, mf, type = "beg") # beg time
  standata$t_end <- make_t_for_stan(formula, mf, type = "end") # end time
  standata$t_gap <- make_t_for_stan(formula, mf, type = "gap") # gap time

  # event indicator for each row of data
  standata$d <- make_d_for_stan(formula, mf)

  # uncensored event times
  tt <- standata$t_end[standata$d == 1]

  # knot locations for fpm baseline hazard
  knots <- get_knots(tt, df = df, iknots = iknots, bknots = bknots)
  iknots <- knots$iknots # internal knot locations
  bknots <- knots$bknots # boundary knot locations
  df     <- knots$df     # validated degrees of freedom

  # basis terms for fpm baseline hazard
  standata$fpm_x_beg <- rcs(standata$t_beg, iknots = iknots, bknots = bknots)
  standata$fpm_x_end <- rcs(standata$t_end, iknots = iknots, bknots = bknots)

  # first derivative of basis terms for fpm baseline hazard
  standata$fpm_dx_beg <- drcs(standata$t_beg, iknots = iknots, bknots = bknots)
  standata$fpm_dx_end <- drcs(standata$t_end, iknots = iknots, bknots = bknots)

  # design matrices for linear predictor
  x <- make_x_for_stan(formula, mf) # fe predictor matrix
  standata$nrows <- x$N
  standata$K     <- x$K
  standata$x     <- x$x
  standata$xbar  <- x$xbar
  #standata$z <- make_z_for_stan(formula, mf) # re predictor matrix
  #standata$g <- make_g_for_stan(formula, mf) # re group ids (for each row)

  #-----------
  # Fit model
  #-----------

  stanfit  <- stanmodels$surv
  stanpars <- pars_to_monitor(standata)
  if (algorithm == "sampling") {
    args <- set_sampling_args(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      prior  = NULL,
      user_dots = dots,
      user_adapt_delta = adapt_delta,
      user_max_treedepth = max_treedepth,
      init = init,
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
  } else {
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm # meanfield or fullrank vb
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }
  check_stanfit(stanfit)

  fit <- nlist(stanfit, formula, family, weights, M, cnms, flevels, n_grps, n_yobs,
               algorithm, terms, glmod = y_mod, data, prior.info = prior_info,
               stan_function = "stan_haz", call = match.call(expand.dots = TRUE))

  out <- stansurv(fit)
  return(out)
}
