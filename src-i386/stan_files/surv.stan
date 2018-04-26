functions {

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
  * @param prior_mu,prior_sd Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  vector make_beta(vector z_beta, int prior_dist, vector prior_mu,
                   vector prior_sd, vector prior_df, real global_prior_sd,
                   real[] global, vector[] local, real[] ool, vector[] mix,
                   real[] aux, int family, real slab_sd, real[] caux) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_sd + prior_mu;
/*    else if (prior_dist == 2) for (k in 1:rows(prior_mu)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_sd[k] + prior_mu[k];
    }
    else if (prior_dist == 3) {
      real c2 = square(slab_sd) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_sd, aux[1], c2);
      else
        beta = hs_prior(z_beta, global, local, global_prior_sd, 1, c2);
    }
    else if (prior_dist == 4) {
      real c2 = square(slab_sd) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_sd, aux[1], c2);
      else
        beta = hsplus_prior(z_beta, global, local, global_prior_sd, 1, c2);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mu + prior_sd .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mu + ool[1] * prior_sd .* sqrt(2 * mix[1]) .* z_beta;
*/    return beta;
  }

  /**
  * Log-prior for coefficients
  *
  * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_sd Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_sd,
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
      // unorthodox useage of prior_sd as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_sd, 0.5 * prior_sd);
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
  matrix[nrows,df] fpm_x_beg;  // design matrix (fpm basis terms)
  matrix[nrows,df] fpm_x_end;  // design matrix (fpm basis terms)
  matrix[nrows,df] fpm_dx_beg; // design matrix (deriv of fpm basis terms)
  matrix[nrows,df] fpm_dx_end; // design matrix (deriv of fpm basis terms)

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
  real                prior_mu_for_exp_scale; // exponential model
  real<lower=0>       prior_sd_for_exp_scale;
  real<lower=0>       prior_df_for_exp_scale;
  real                prior_mu_for_wei_shape; // weibull model
  real<lower=0>       prior_sd_for_wei_shape;
  real<lower=0>       prior_df_for_wei_shape;
  real                prior_mu_for_wei_scale;
  real<lower=0>       prior_sd_for_wei_scale;
  real<lower=0>       prior_df_for_wei_scale;
  vector[df]          prior_mu_for_fpm_coefs; // fpm model
  vector<lower=0>[df] prior_sd_for_fpm_coefs;
  vector<lower=0>[df] prior_df_for_fpm_coefs;
}

transformed data {
  int<lower=0> hs = get_nvars_for_hs(prior_dist);
}

parameters {

  // primitive log hazard ratios
  vector[K] z_beta;

  // unscaled basehaz parameters
  real<lower=0> z_exp_scale[dist == 1]; // unscaled scale par (exponential model)
  real<lower=0> z_wei_shape[dist == 2]; // unscaled shape par (weibull model)
  real<lower=0> z_wei_scale[dist == 2]; // unscaled scale par (weibull model)
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
  beta = make_beta(z_beta, prior_dist, prior_mu,
                   prior_sd, prior_df, global_prior_sd,
                   global, local, ool, mix, rep_array(1.0, 0), 0,
                   slab_sd, caux);

  // define basehaz parameters
  if (dist == 1) { // exponential model
    exp_scale[1] = z_exp_scale[1] * prior_sd_for_exp_scale;
  }
  else if (dist == 2) { // weibull model
    wei_shape[1] = z_wei_shape[1] * prior_sd_for_wei_shape;
    wei_scale[1] = z_wei_scale[1] * prior_sd_for_wei_scale;
  }
  else if (dist == 3) { // fpm model
    fpm_coefs[1] = prior_mu_for_fpm_coefs +
      z_fpm_coefs[1] .* prior_sd_for_fpm_coefs;
  }
}

model {

	vector[nrows] eta; // linear predictor
  vector[nrows] log_basehaz; // log basehaz at t_end for each row of data
  vector[nrows] log_basesurv_beg; // log basesurv at t_beg for each row of data
  vector[nrows] log_basesurv_end; // log basesurv at t_end for each row of data
  vector[nrows] log_haz;  // log haz at t_end for each row of data
  vector[nrows] log_surv; // log surv at t_end for each row of data

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
		for (n in 1:nrows) {
			log_basesurv_beg[n] = - wei_scale[1] * (t_beg[n] ^ wei_shape[1]);
			log_basesurv_end[n] = - wei_scale[1] * (t_end[n] ^ wei_shape[1]);
		}
	}
	else if (dist == 3) { // fpm model
	  log_basehaz = - log(t_end) + log(fpm_dx_end * fpm_coefs[1]) + (fpm_x_end * fpm_coefs[1]);
		log_basesurv_beg = - exp(fpm_x_beg * fpm_coefs[1]);
		log_basesurv_end = - exp(fpm_x_end * fpm_coefs[1]);
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
