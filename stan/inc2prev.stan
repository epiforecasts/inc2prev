functions {
#include gaussian_process.stan
#include rt.stan
#include convolve.stan
#include observed_in_window.stan
#include ab.stan
#include generated_quantities.stan
#include rng.stan
}

data {
  int<lower = 0> ut; // initial period (before data starts)
  int<lower = 0> t; // number of time points to model
  int<lower = 0> n; // number of time series to model
  int<lower = 0> n_ab; // number of antibody time series to model
  array[n_ab] int<lower = 0> ab_index; // antibody dynamics indices
  int<lower = 0> obs; // number of prevalence observations
  int<lower = 0> ab_obs; // presence/absence of antibody observations (1 = present; 0 = absent)
  array[n] vector<lower = 0, upper = 1>[obs] prev; // observed positivity prevalence
  array[n] vector<lower = 0>[obs] prev_sd2; // squared standard deviation of observed positivity prevalence
  array[n_ab] vector<lower = 0, upper = 1>[ab_obs] ab; // observed antibody posivitiy prevalence
  array[n_ab] vector<lower = 0>[ab_obs] ab_sd2; // squared standard deviation of observed antibody prevalence
  array[obs] int<lower = 0> prev_stime; // starting times of positivity prevalence observations
  array[obs] int <lower = 0>prev_etime; // end times of positivity prevalence observations
  array[ab_obs] int<lower = 0> ab_stime; // starting times of antibody prevalence observations
  array[ab_obs] int<lower = 0> ab_etime; // end times of antibody prevalence observations
  array[n] vector<lower = 0>[n_ab > 0 ? t : 0] vacc; // vaccinations
  int<lower = 0> pbt; // maximum detection time
  vector<lower = 0, upper = 1>[pbt] prob_detect_mean; // at each time since infection, probability of detection
  vector<lower = 0, upper = 1>[pbt] prob_detect_sd; // at each time since infection, tandard deviation of probability of detection
  real<lower = 0> lengthscale_alpha; // alpha for gp lengthscale prior
  real<lower = 0> lengthscale_beta;  // beta for gp lengthscale prior
  int<lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  int<lower = 0, upper = 2> diff_order; // Order of differencing to use for the GP (0 = cases -> mean, 1 = growth_t -> 0, 2 = growth_t -> growth_{t-1})
  array[2] real<lower = 0> gtm; // mean and standard deviation (sd) of the mean generation time
  array[2] real<lower = 0> gtsd; // mean and sd of the sd of the generation time
  int<lower = 0> gtmax; // maximum number of days to consider for the generation time
  array[n] real init_inc_mean; // Mean initial/mean incidence (logit)
  array[n_ab > 0 ? n : 0] real<lower = 0> init_ab_mean; // mean estimate of initial antibody prevalence
  array[n_ab > 0 ? n : 0] real<lower = 0> init_ab_sd; // sd of estimate of initial antibody prevalence
  array[2] real pbeta; // Mean and sd for prior proportion that seroconvert
  array[2] real pgamma_mean; // Means for prior infection and vaccine waning
  array[2] real pgamma_sd; // Sds for prior infection and vaccine waning
  array[2] real pdelta; // Mean and sd for prior vaccine efficacy
  int linf_ab_delay; // Length of the delay between infection and ab
  // PMF of the delay between infection and ab
  vector<lower = 0>[linf_ab_delay] inf_ab_delay; 
  int<lower = 0> lvacc_ab_delay; // Length of the delay between vaccination and ab
  // PMF of the  delay between vaccination and ab
  vector<lower = 0>[lvacc_ab_delay] vacc_ab_delay; 
  int<lower = 0, upper = 1> prev_likelihood; // Should the likelihood for prevalence data be included
  int<lower = 0, upper = 1> ab_likelihood; // Should the likelihood for antibody data be included
}

transformed data {
  array[n] vector<lower = 0, upper = 1>[n_ab > 0 ? t : 0] vacc_with_ab;
  // set up approximate gaussian process
  matrix[t - diff_order, M] PHI = setup_gp(M, L, t - diff_order);
  // Calculate vaccinations with the potential to have antibodies
  if (n_ab > 0) {
    for (i in 1:n) {
      vacc_with_ab[i] = convolve(vacc[i], vacc_ab_delay);
    }
  }
}

parameters {
  array[n] real<lower = 0> rho; // length scale of gp
  array[n] real<lower = 0> alpha; // scale of gp
  array[n] vector[M] eta; // eta of gp
  array[n] real init_inc; // Initial infections
  array[n] vector[diff_order] init_growth;
  real<lower = 0> sigma; // observation error
  array[n_ab > 0 ? 1 : 0] real<lower = 0> ab_sigma; // observation error
  vector<lower = 0, upper = 1>[pbt] prob_detect; // probability of detection as a function of time since infection
  vector[n_ab > 0 ? 1 : 0] logit_beta; // proportion that seroconvert
  vector[n_ab > 0 ? 2 : 0] logit_gamma; // antibody waning (inf & vac)
  vector[n_ab > 0 ? 1 : 0] logit_delta; // vaccine efficacy
  vector<lower = 0, upper = 1>[n_ab > 0 ? n : 0] init_dab; // initial proportion with antibodies
  vector[n_ab > 0 ? 1 : 0] k; // Potential loss of efficacy from new infections in already seropositive people
  vector[n_ab > 0 ? 1 : 0] l; // Potential loss of efficacy from new doses being administered to already seropositive people
}

transformed parameters {
  array[n] vector[t] gp; // value of gp at time t + initialisation
  array[n_ab] vector[t] infs_with_potential_abs; // Infections with the potential to have ab
  array[n] vector<lower = 0, upper = 1>[t] infections; // incident infections at time t
  array[n] vector<lower = 0>[t] dcases; // detectable cases at time t
  array[n_ab] vector<lower = 0>[t] dab; // proportion of individuals with antibodies at time t
  array[n] vector<lower = 0>[obs] odcases;
  array[n_ab] vector<lower = 0>[ab_obs] odab;
  array[n] vector[obs] combined_sigma;
  array[n_ab] vector[ab_obs] combined_ab_sigma;
  vector[n_ab > 0 ? 1 : 0] beta = inv_logit(logit_beta); // tranformation to natural scale
  vector[n_ab > 0 ? 2 : 0] gamma = inv_logit(logit_gamma); // tranformation to natural scale
  vector[n_ab > 0 ? 1 : 0] delta = inv_logit(logit_delta); // tranformation to natural scale
  

  // update gaussian process
  for (i in 1:n) {
    gp[i][(1 + diff_order):t] = update_gp(PHI, M, L, alpha[i], rho[i], eta[i], 0);
    // setup differencing of the GP
    if (diff_order) {
      gp[i][1:diff_order] = init_growth[i];
      for (j in 1:diff_order) {
        gp[i] = cumulative_sum(gp[i]);
      }
    }

    // relative probability of infection
    // inc_init is the mean incidence
    infections[i] = inv_logit(init_inc[i] + gp[i]);

    // calculate detectable cases
    dcases[i] = convolve(infections[i], prob_detect);
    // calculate observed detectable cases
    odcases[i] = observed_in_window(dcases[i], prev_stime, prev_etime, ut, obs);
    //combined standard error
    combined_sigma[i] = sqrt(square(sigma) + prev_sd2[i]);
  }

  //calculate infections with potential to have antibodies
  for (i in 1:n_ab) {
    infs_with_potential_abs[i] = convolve(infections[ab_index[i]], inf_ab_delay);
    // calculate detectable antibodies
    dab[i] = detectable_antibodies(infs_with_potential_abs[i],
                                   vacc_with_ab[ab_index[i]], 
                                   beta[1], gamma, delta[1], k[1], l[1],
                                   init_dab[ab_index[i]], t);
    // calculate observed detectable antibodies
    odab[i] = observed_in_window(dab[i], ab_stime, ab_etime, ut, ab_obs);
    //combined standard error
    combined_ab_sigma[i] = sqrt(square(ab_sigma[1]) + ab_sd2[i]);
  }
}

model {
  // gaussian process priors
  for (i in 1:n) {
    rho[i] ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
    alpha[i] ~ normal(0, pow(0.1, diff_order)) T[0,];
    for (j in 1:M) {
      eta[i, j] ~ std_normal();
    }
    init_inc[i] ~ normal(init_inc_mean[i], 2);
  }

  // Initial infections
  for (i in 1:n) {
    for (j in 1:diff_order) {
      init_growth[i, j] ~ normal(0, 0.25);
    }
  }

  // prevalence observation model
  for (i in 1:pbt) {
    prob_detect[i] ~ normal(prob_detect_mean[i], prob_detect_sd[i]) T[0, 1];
  }

  // Priors for antibody model
  if (n_ab > 0) {
    init_dab ~ normal(init_ab_mean, init_ab_sd);
    logit_beta ~ normal(pbeta[1], pbeta[2]);
    logit_gamma ~ normal(pgamma_mean, pgamma_sd); 
    logit_delta ~ normal(pdelta[1], pdelta[2]);
    k ~ lognormal(0, 0.1);
    l ~ lognormal(0, 0.1);
    ab_sigma[1] ~ normal(0.025, 0.025) T[0,];
  }

  sigma ~ normal(0.005, 0.0025) T[0,];
 
  if (prev_likelihood) {
    for (i in 1:n) {
      prev[i] ~ normal(odcases[i], combined_sigma[i]);
    }
  }
  if (ab_likelihood && ab_obs) {
    for (i in 1:n_ab) {
      ab[i] ~ normal(odab[i], combined_ab_sigma[i]);
    }
  }
}

generated quantities {
  array[n] vector<lower = 0>[t - ut] R;
  array[n] vector[t - ut] r;
  array[n, obs] real<lower = 0, upper = 1> est_prev;
  array[n_ab > 0 ? n : 0, ab_obs] real<lower = 0, upper = 1> est_ab;
  array[n_ab > 0 ? n : 0] vector<lower = 0, upper = 1>[t] gen_dab;
  array[n_ab > 0 ? n : 0] vector<lower = 0, upper = 1>[t] gen_infs_with_potential_abs; // Infections with the potential to have ab
  array[n_ab > 0 ? n : 0] vector<lower = 0, upper = 1>[ab_obs] gen_odab;
  real<lower = 0> gtm_sample;
  real<lower = 0> gtsd_sample;

  {
    // sample generation time
    gtm_sample = normal_rng(gtm[1], gtm[2]);
    gtsd_sample = normal_rng(gtsd[1], gtsd[2]);

    // get cumulative incidence
    for (i in 1:n) {
      // sample estimated prevalence
      for (j in 1:obs) {
        est_prev[i, j] = normal_lub_rng(odcases[i, j], sigma, 0, 1);
      }
      // calculate Rt using infections and generation time
      R[i] = calculate_Rt(infections[i], ut, gtm_sample, gtsd_sample, gtmax, 1);
      // calculate growth
      r[i] = calculate_growth(infections[i], ut);
    }

    if (n_ab > 0) {
      for (i in 1:n) {
        gen_infs_with_potential_abs[i] = convolve(infections[i], inf_ab_delay);
        // calculate detectable antibodies
        gen_dab[i] = detectable_antibodies(gen_infs_with_potential_abs[i], 
                                           vacc_with_ab[i], beta[1],
                                           gamma, delta[1], k[1], l[1], init_dab[i], t);
        // calculate observed detectable antibodies
        gen_odab[i] = observed_in_window(gen_dab[i], ab_stime, ab_etime, ut, ab_obs);
        // estimate with error
        for (j in 1:ab_obs) {
          est_ab[i, j] = normal_lub_rng(gen_odab[i, j], ab_sigma[1], 0, 1);
	}
      }
    }
  }
}
