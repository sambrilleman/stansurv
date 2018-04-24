functions {

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
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
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
    return beta;
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
    else if (prior_dist == 3) { // hs
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
    /* else prior_dist is 0 and nothing is added */
  }

}

data {

  // data for event model
  int<lower=0> K;             // num. cols in predictor matrix
  int<lower=0> nrows;         // num. rows in predictor matrix
  int<lower=0> npats;         // num. individuals
  int<lower=0> nevents;       // num. events (ie. not censored)
  int<lower=0> df;            // df for baseline hazard
  int<lower=0> df_tde[K];     // df for time-dependent hazard ratios
  vector[nrows] t_beg;        // beg time for each row of data
  vector[nrows] t_end;        // end time for each row of data
  vector[nrows] t_gap;        // gap time for each row of data
  vector[nrows] d;            // event indicator for each row of data
  matrix[K > 0 ? nrows : 0, K] x; // predictor matrix
  matrix[nrows,df] fpm_x_beg; // design matrix (basis terms) for basehaz
  matrix[nrows,df] fpm_x_end; // design matrix (basis terms) for basehaz

  // survival distribution:
  //   1 = exponential
  //   2 = weibull
  //   3 = fpm
  int<lower=1,upper=3> dist;

	// flag to draw from prior predictive distribution
	int<lower=0,upper=1> prior_PD; // 0 = no, 1 = yes

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
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_exp_scale;
  int<lower=0,upper=3> prior_dist_for_wei_shape;
  int<lower=0,upper=3> prior_dist_for_wei_scale;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> prior_dist_for_fpm_coefs;

  // hyperparameter (log hazard ratios), set to 0 if there is no prior
  vector[K]           prior_mu;
  vector<lower=0>[K]  prior_sd;
  vector<lower=0>[K]  prior_df;
  real<lower=0>       global_prior_sd; // for hs priors only
  real<lower=0>       global_prior_df;
  real<lower=0>       slab_sd;
  real<lower=0>       slab_df;

  // hyperparameters (basehaz pars), set to 0 if there is no prior
  real            prior_mu_for_exp_scale; // exponential model
  real<lower=0>   prior_sd_for_exp_scale;
  real<lower=0>   prior_df_for_exp_scale;
  real            prior_mu_for_wei_shape; // weibull model
  real<lower=0>   prior_sd_for_wei_shape;
  real<lower=0>   prior_df_for_wei_shape;
  real            prior_mu_for_wei_scale;
  real<lower=0>   prior_sd_for_wei_scale;
  real<lower=0>   prior_df_for_wei_scale;
  vector          prior_mu_for_fpm_coefs; // fpm model
  vector<lower=0> prior_sd_for_fpm_coefs;
  vector<lower=0> prior_df_for_fpm_coefs;
}

parameters {

  // primitive log hazard ratios
  vector[K] z_beta;

  // unscaled basehaz parameters
  real<lower=0> z_exp_scale[dist == 1]; // unscaled scale par (exponential model)
  real<lower=0> z_wei_shape[dist == 2]; // unscaled shape par (weibull model)
  real<lower=0> z_wei_lscale[dist == 2]; // unscaled scale par (weibull model)
  vector[df]    z_fpm_coefs[dist == 3]; // unscaled spline coefs (fpm model)

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
  real<lower=0> exp_scale[dist == 1]; // scale par (exp dist)
  real<lower=0> wei_shape[dist == 2]; // shape par (weibull dist)
  real<lower=0> wei_scale[dist == 2]; // scale par (weibull dist)
  vector[df]    fpm_coefs[dist == 3]; // spline coefs (fpm model)

  // define log hazard ratios
  beta = make_beta(z_beta, prior_dist, prior_mean,
                   prior_scale, prior_df, lobal_prior_scale,
                   global, local, ool, mix, rep_array(1.0, 0), 0,
                   slab_scale, caux);

  // define basehaz parameters
  if (dist == 1) { // exponential model
    exp_scale = z_exp_scale * prior_scale_for_exp_scale;
  }
  else if (dist == 2) { // weibull model
    wei_shape = z_wei_shape * prior_scale_for_wei_shape;
    wei_scale = z_wei_scale * prior_scale_for_wei_scale;
  }
  else if (dist == 3) { // fpm model
    fpm_coefs = prior_mu_for_fpm_coefs +
      z_fpm_coefs .* prior_scale_for_fpm_coefs;
  }
}

model {

	vector[nrows] eta; // linear predictor
  vector[nrows] log_basehaz; // log basehaz at t_end for each row of data
  vector[nrows] log_basesurv_beg; // log survival at t_beg for each row of data
  vector[nrows] log_basesurv_end; // log survival at t_end for each row of data

	// linear predictor
	if (K > 0) {
	  eta = x * beta;
	}
	else {
	  eta = rep_vector(0.0, nrows);
	}

  // log basehaz and log basesurv for each row of data
  if (dist == 1) { // exponential model
    log_basehaz = rep_vector(log(exp_scale[1]), nrows);
		log_basesurv_beg = - exp_scale[1] * t_beg;
		log_basesurv_end = - exp_scale[1] * t_end;
	}
	else if (dist == 2) { // weibull model
	  log_basehaz = log(wei_shape[1]) + log(wei_scale[1]) + t_end * (wei_shape[1] - 1);
		log_basesurv_beg = - wei_scale[1] * (t_beg ^ wei_shape[1]);
		log_basesurv_end = - wei_scale[1] * (t_end ^ wei_shape[1]);
	}
	else if (dist == 3) { // fpm model
	  log_basehaz = - log(t_end) + log(fpm_dx_end * fpm_coefs) + (fpm_x_end * fpm_coefs);
		log_basesurv_beg = - exp(fpm_x_beg * fpm_coefs);
		log_basesurv_end = - exp(fpm_x_end * fpm_coefs);
	}

  // log hazard for each row of data
  log_haz = log_basehaz + eta;

	// log survival for each row of data (allows for delayed entry)
	log_surv = (log_basesurv_end - log_basesurv_beg) .* exp(eta);

  // log likelihood for event model
  if (prior_PD == 0) { // unweighted log likelihood
    target += sum(d .* log_haz) + sum(log_surv);
  }

	// log priors for coefficients
	beta_lp(z_beta, prior_dist, prior_sd, prior_df, global_prior_df,
	        local, global, mix, ool, slab_df, caux);

	// log priors for baseline hazard parameters

}
