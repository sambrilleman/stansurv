#' Bayesian proportional hazards regression
#'
#' Bayesian inference for proportional hazards regression models. The user
#' can specify a variety of standard parametric distributions for the
#' baseline hazard, or a Royston-Parmar flexible parametric model.
#'
#' @export
#'
#' @examples
#' pbc2 <- survival::pbc
#' pbc2 <- pbc2[!is.na(pbc2$trt),]
#' pbc2$status <- as.integer(pbc2$status > 0)
#' m1 <- stan_surv(survival::Surv(time, status) ~ trt, data = pbc2)
#'
#' df <- flexsurv::bc
#' m2 <- stan_surv(survival::Surv(rectime, censrec) ~ group,
#'                 data = df, cores = 1, chains = 1, iter = 2000,
#'                 basehaz = "fpm", iknots = c(6.594869,  7.285963 ),
#'                 degree = 2, prior_aux = normal(0, 2, autoscale = F))
#'
stan_surv <- function(formula, data, basehaz = "fpm", timescale = "log",
                      df = 5L, degree = 3L, iknots = NULL, bknots = NULL,
                      prior = normal(), prior_intercept = normal(),
                      prior_aux = list(), prior_PD = FALSE,
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

  #----------------
  # Construct data
  #----------------

  standata <- list()

  # model data frame
  mf <- data

  #----- dimensions, response, predictor matrix

  # time variable for each row of data
  standata$t_beg <- make_t(formula, mf, type = "beg") # beg time
  standata$t_end <- make_t(formula, mf, type = "end") # end time
  standata$t_gap <- make_t(formula, mf, type = "gap") # gap time

  # event indicator for each row of data
  standata$d <- make_d(formula, mf)

  # design matrices for linear predictor
  x <- make_x(formula, mf) # fe predictor matrix
  standata$nrows <- x$N
  standata$K     <- x$K
  standata$x     <- x$x
  standata$xbar  <- x$xbar
  #standata$z <- make_z(formula, mf) # re predictor matrix
  #standata$g <- make_g(formula, mf) # re group ids (for each row)

  #----- time-dependent effects (i.e. non-proportional hazards)

  # degrees of freedom for time-dependent effects
  standata$df_tde <- aa(rep(0L, standata$K)) # not implemented yet

  #----- baseline hazard

  ok_basehaz <- c("exponential", "weibull", "fpm", "fpm2")
  basehaz <- handle_basehaz(basehaz, df = df, degree = degree,
                            iknots = iknots, bknots = bknots,
                            t_beg = standata$t_beg, t_end = standata$t_end,
                            d = standata$d, ok_basehaz = ok_basehaz,
                            timescale = timescale)
  standata$type <- ai(basehaz$type)
  standata$df   <- ai(basehaz$df)
  #standata$basehaz_x_beg  <- make_basehaz_x(standata$t_beg, basehaz = basehaz)
  #standata$basehaz_dx_beg <- make_basehaz_x(standata$t_beg, basehaz = basehaz, deriv = TRUE)
  standata$basehaz_x_beg  <- matrix(0, standata$nrows, standata$df)
  standata$basehaz_dx_beg <- matrix(0, standata$nrows, standata$df)
  standata$basehaz_x_end  <- make_basehaz_x(standata$t_end, basehaz = basehaz,
                                            timescale = timescale)
  standata$basehaz_dx_end <- make_basehaz_x(standata$t_end, basehaz = basehaz,
                                            timescale = timescale, deriv = TRUE)
  standata$has_intercept  <- ai(has_intercept(basehaz))

  #----- priors and hyperparameters

  # priors
  user_prior_stuff <- prior_stuff <-
    handle_prior(prior, nvars = x$K,
                 default_scale = 2,
                 ok_dists = ok_dists())

  user_prior_intercept_stuff <- prior_intercept_stuff <-
    handle_prior(prior_intercept, nvars = 1,
                 default_scale = 2,
                 ok_dists = ok_dists_for_intercept())

  user_prior_aux_stuff <- prior_aux_stuff <-
    handle_prior(prior_aux, nvars = basehaz$df,
                 default_scale = get_default_aux_scale(basehaz$name),
                 ok_dists = ok_dists_for_aux())

  # autoscaling of priors
  prior_stuff           <- autoscale_prior(prior_stuff, predictors = x$x)
  prior_intercept_stuff <- autoscale_prior(prior_intercept_stuff)
  prior_aux_stuff       <- autoscale_prior(prior_aux_stuff)

  # priors
  standata$prior_dist              <- prior_stuff$prior_dist
  standata$prior_dist_for_intercept<- prior_intercept_stuff$prior_dist
  standata$prior_dist_for_aux      <- prior_aux_stuff$prior_dist

  # hyperparameters
  standata$prior_mean               <- prior_stuff$prior_mean
  standata$prior_scale              <- prior_stuff$prior_scale
  standata$prior_df                 <- prior_stuff$prior_df
  standata$prior_mean_for_intercept <- c(prior_intercept_stuff$prior_mean)
  standata$prior_scale_for_intercept<- c(prior_intercept_stuff$prior_scale)
  standata$prior_df_for_intercept   <- c(prior_intercept_stuff$prior_df)
  standata$prior_scale_for_aux      <- prior_aux_stuff$prior_scale
  standata$prior_df_for_aux         <- prior_aux_stuff$prior_df
  standata$global_prior_scale       <- prior_stuff$global_prior_scale
  standata$global_prior_df          <- prior_stuff$global_prior_df
  standata$slab_df                  <- prior_stuff$slab_df
  standata$slab_scale               <- prior_stuff$slab_scale

  #----- additional flags

  standata$prior_PD <- ai(prior_PD)
  standata$delayed <- ai(!all_zero(standata$t_beg))
  standata$npats <- standata$nevents <- 0L # not currently used

  #-----------
  # Fit model
  #-----------

  stanfit  <- stanmodels$surv
  stanpars <- pars_to_monitor(standata)
  if (algorithm == "sampling") {
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      show_messages = FALSE
    )
#    args <- set_sampling_args(
#      object = stanfit,
#      data   = standata,
#      pars   = stanpars,
#      prior  = NULL,
#      user_dots = dots,
#      user_adapt_delta = adapt_delta,
#      user_max_treedepth = max_treedepth,
#      init = init,
#      show_messages = FALSE)
    args[names(dots)] <- dots
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
  #check_stanfit(stanfit)

  fit <- nlist(stanfit, formula, data, basehaz, algorithm,
               stan_function = "stan_surv", call = match.call(expand.dots = TRUE))

  #out <- stansurv(fit)
  return(fit)
}
