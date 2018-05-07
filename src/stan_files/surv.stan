functions {

  /**
  * Return the lower bound for the baseline hazard parameters
  *
  * @param type An integer indicating the type of baseline hazard
  * @return A real
  */
  real basehaz_coefs_lb(int type) {
	  real lb = (type == 4) ? negative_infinity() : 0;
    return lb;
	}
	
  /**
  * Return the required number of local hs parameters
  *
  * @param prior_dist An integer indicating the prior distribution
  * @return An integer
  */
  int get_nvars_for_hs(int prior_dist) {
    int hs = 0;
    if (prior_dist == 3) hs = 2;
    else if (prior_dist == 4) hs = 4;
    return hs;
  }

  /**
  * Scale the primitive population level parameters based on prior information
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  vector make_beta(vector z_beta, int prior_dist, vector prior_mean,
                   vector prior_scale, vector prior_df, real global_prior_scale,
                   real[] global, vector[] local, real[] ool, vector[] mix,
                   real[] aux, int family, real slab_scale, real[] caux) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
/*    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 4) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + ool[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
*/    return beta;
  }

  /**
  * Log-prior for coefficients
  *
  * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_scale Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df, real global_prior_df, vector[] local,
               real[] global, vector[] mix, real[] one_over_lambda,
               real slab_df, real[] caux) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
/*    else if (prior_dist == 3) { // hs
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 4) { // hs+
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 5) { // laplace
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
      target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      target += normal_lpdf(z_beta | 0, 1);
    }
*/    /* else prior_dist is 0 and nothing is added */
  }

  /**
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector, the unscaled auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param df Vector, the df for the prior distribution
  * @return nothing
  */
  void aux_lp(vector aux_unscaled, int dist, vector df) {
    if (dist > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
  }

}

data {

  // dimensions, response, predictor matrix
  int<lower=0> K;          // num. cols in predictor matrix
  int<lower=0> nrows;      // num. rows in predictor matrix
  int<lower=0> npats;      // num. individuals
  int<lower=0> nevents;    // num. events (ie. not censored)
  int<lower=0> df;         // df for baseline hazard
  int<lower=0> df_tde[K];  // df for time-dependent hazard ratios
  vector[nrows] t_beg;     // beg time for each row of data
  vector[nrows] t_end;     // end time for each row of data
  vector[nrows] d;         // event indicator for each row of data
  matrix[nrows,K] x;       // predictor matrix
  int<lower=0,upper=1> delayed; // flag for delayed entry
  int<lower=0,upper=1> has_intercept; // flag for intercept

  // baseline hazard type:
  //   1 = exponential
  //   2 = weibull
  //   3 = fpm  (fpm on cum haz scale)
	//   4 = fpm2 (fpm on log cum haz scale)
  int<lower=1,upper=4> type;

  // data for baseline hazard
  //   exp model:
  //     basehaz_x:  not used
  //     basehaz_dx: not used
  //   wei model:
  //     basehaz_x:  log times
  //     basehaz_dx: not used
  //   fpm model:
  //     basehaz_x:  I-spline basis terms
  //     basehaz_dx: I-spline derivative of basis terms
  matrix[nrows,df] basehaz_x_beg;
  matrix[nrows,df] basehaz_x_end;
  matrix[nrows,df] basehaz_dx_beg;
  matrix[nrows,df] basehaz_dx_end;

  // flag to draw from prior predictive distribution
  int<lower=0,upper=1> prior_PD; // 1 = yes

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> prior_dist;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> prior_dist_for_intercept;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;

  // hyperparameter (log hazard ratios), set to 0 if there is no prior
  vector[K]           prior_mean;
  vector<lower=0>[K]  prior_scale;
  vector<lower=0>[K]  prior_df;
  real<lower=0>       global_prior_scale; // for hs priors only
  real<lower=0>       global_prior_df;
  real<lower=0>       slab_scale;
  real<lower=0>       slab_df;

  // hyperparameters (intercept), set to 0 if there is no prior
  real                prior_mean_for_intercept;
  real<lower=0>       prior_scale_for_intercept;
  real<lower=0>       prior_df_for_intercept;

  // hyperparameters (basehaz pars), set to 0 if there is no prior
  vector<lower=0>[df] prior_scale_for_aux;
  vector<lower=0>[df] prior_df_for_aux;
}

transformed data {
  int<lower=0> hs = get_nvars_for_hs(prior_dist);
}

parameters {

  // primitive log hazard ratios
  vector[K] z_beta;

  // intercept
  real gamma[has_intercept == 1];

  // unscaled basehaz parameters
  //   exp model: df = 0, ie. no aux parameter
  //   wei model: df = 1, ie. 1 shape parameter
  //   fpm model: df = number of basis terms, ie. I-spline coefs
  vector<lower=basehaz_coefs_lb(type)>[df] z_basehaz_coefs;

  // parameters for priors
  real<lower=0> global[hs];
  vector<lower=0>[hs > 0 ? K : 0] local[hs];
  real<lower=0> caux[hs > 0];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> ool[prior_dist == 6];
}

transformed parameters {

  // log hazard ratios
  vector[K] beta;

  // basehaz parameters
  vector[df] basehaz_coefs;

  // define log hazard ratios
  beta = make_beta(z_beta, prior_dist, prior_mean,
                   prior_scale, prior_df, global_prior_scale,
                   global, local, ool, mix, rep_array(1.0, 0), 0,
                   slab_scale, caux);

  // define basehaz parameters
  if (type > 1) {
    basehaz_coefs = z_basehaz_coefs .* prior_scale_for_aux;
  }
}

model {

  vector[nrows] eta;      // linear predictor
  vector[nrows] log_haz;  // log haz  at t_end for each row of data
  vector[nrows] log_surv_end; // log surv at t_end for each row of data
  vector[nrows] log_surv_beg; // log surv at t_beg for each row of data

  // linear predictor
  if (K > 0) {
    eta = x * beta;
  }
  else {
    eta = rep_vector(0.0, nrows);
  }

  // add intercept
  if (has_intercept == 1) {
    eta = eta + gamma[1];
  }

  // log basehaz and log basesurv for each row of data
  if (type == 1) { // exponential model
    vector[nrows] exp_eta = exp(eta);
    log_haz = eta;
    log_surv_end = - t_end .* exp_eta;
    if (delayed == 1)
      log_surv_beg = - t_beg .* exp_eta;
  }
  else if (type == 2) { // weibull model
    vector[nrows] exp_eta = exp(eta);
    log_haz = log(basehaz_coefs[1]) + basehaz_x_end[,1] * (basehaz_coefs[1] - 1) + eta;
    for (n in 1:nrows) {
      log_surv_end[n] = - (t_end[n] ^ basehaz_coefs[1]) * exp_eta[n];
      if (delayed == 1)
        log_surv_beg[n] = - (t_beg[n] ^ basehaz_coefs[1]) * exp_eta[n];
    }
  }
  else if (type == 3) { // fpm model -- cum haz scale
    vector[nrows] exp_eta = exp(eta);
    log_haz = log(basehaz_dx_end * basehaz_coefs) + eta;
    log_surv_end = - (basehaz_x_end * basehaz_coefs) .* exp_eta;
    if (delayed == 1)
      log_surv_beg = - (basehaz_x_beg * basehaz_coefs) .* exp_eta;
  }	
  else if (type == 4) { // fpm2 model -- log cum haz scale
    log_haz = - log(t_end) + log(basehaz_dx_end * basehaz_coefs) +
      (basehaz_x_end * basehaz_coefs);
    log_surv_end = - exp(basehaz_x_end * basehaz_coefs + eta);
    if (delayed == 1)
      log_surv_beg = - exp(basehaz_x_beg * basehaz_coefs + eta);
  }
	
  // correct log likelihood for delayed entry
  if (delayed == 1) {
    log_surv_end = log_surv_end - log_surv_beg;
  }

  // log likelihood for event model
  if (prior_PD == 0) { // unweighted log likelihood
    target += d .* log_haz + log_surv_end;
  }

  // log priors for coefficients
  beta_lp(z_beta, prior_dist, prior_scale, prior_df, global_prior_df,
          local, global, mix, ool, slab_df, caux);

  // log prior for intercept
  if (has_intercept == 1) {
    target += normal_lpdf(gamma[1] | prior_mean_for_intercept,
                                     prior_scale_for_intercept);
  }

  // log priors for baseline hazard parameters
  if (type > 1) {
    aux_lp(z_basehaz_coefs, prior_dist_for_aux, prior_df_for_aux);
  }

}
