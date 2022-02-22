functions {
#include gaussian_process.stan
#include rt.stan
#include convolve.stan
#include observed_in_window.stan
#include ab.stan
#include generated_quantities.stan
}

data {
  int ut; // initial period (before data starts)
  int t; // number of time points to model
  int n; // number of time series to model
  int n_ab; // number of antibody time series to model
  int ab_index[n_ab]; // antibody dynamics indices
  int obs; // number of prevalence observations
  int ab_obs; // presence/absence of antibody observations (1 = present; 0 = absent)
  matrix[n, obs] prev; // observed positivity prevalence
  matrix[n, obs] prev_sd2; // squared standard deviation of observed positivity prevalence
  matrix[n_ab, ab_obs] ab; // observed antibody posivitiy prevalence
  matrix[n_ab, ab_obs] ab_sd2; // squared standard deviation of observed antibody prevalence
  int prev_stime[obs]; // starting times of positivity prevalence observations
  int prev_etime[obs]; // end times of positivity prevalence observations
  int ab_stime[ab_obs]; // starting times of antibody prevalence observations
  int ab_etime[ab_obs]; // end times of antibody prevalence observations
  matrix[n, n_ab > 0 ? t : 0] vacc; // vaccinations
  int pbt; // maximum detection time
  vector[pbt] prob_detect_mean; // at each time since infection, probability of detection
  vector[pbt] prob_detect_sd; // at each time since infection, tandard deviation of probability of detection
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  int diff_order; // Order of differencing to use for the GP (0 = cases -> mean, 1 = growth_t -> 0, 2 = growth_t -> growth_{t-1})
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
  real init_inc_mean[n]; // Mean initial/mean incidence (logit)
  real init_ab_mean[n_ab > 0 ? n : 0]; // mean estimate of initial antibody prevalence
  real init_ab_sd[n_ab > 0 ? n : 0]; // sd of estimate of initial antibody prevalence
  real pbeta[2]; // Mean and sd for prior proportion that don't seroconvert
  real pgamma_mean[2]; // Means for prior infection and vaccine waning
  real pgamma_sd[2]; // Sds for prior infection and vaccine waning
  real pdelta[2]; // Mean and sd for prior vaccine efficacy
  int linf_ab_delay; // Length of the delay between infection and ab
  // PMF of the delay between infection and ab
  vector[linf_ab_delay] inf_ab_delay; 
  int lvacc_ab_delay; // Length of the delay between vaccination and ab
  // PMF of the  delay between vaccination and ab
  vector[lvacc_ab_delay] vacc_ab_delay; 
  int prev_likelihood; // Should the likelihood for prevalence data be included
  int ab_likelihood; // Should the likelihood for antibody data be included
}

transformed data {
  matrix[n, n_ab > 0 ? t : 0] vacc_with_ab;
  // set up approximate gaussian process
  matrix[t - diff_order, M] PHI = setup_gp(M, L, t - diff_order);
  // Calculate vaccinations with the potential to have antibodies
  if (n_ab > 0) {
    for (i in 1:n) {
      vacc_with_ab[i] = to_row_vector(convolve(to_vector(vacc[i]), vacc_ab_delay));
    }
  }
}

parameters {
  real <lower = 0> rho[n]; // length scale of gp
  real <lower = 0> alpha[n]; // scale of gp
  matrix[n, M] eta; // eta of gp
  real init_inc[n]; // Initial infections
  matrix[n, diff_order] init_growth;
  real<lower = 0> sigma; // observation error
  real<lower = 0> ab_sigma[n_ab > 0 ? 1 : 0]; // observation error
  vector<lower = 0, upper = 1>[pbt] prob_detect; // probability of detection as a function of time since infection
  vector<lower = 0, upper = 1>[n_ab > 0 ? 1 : 0] beta; // proportion that don't seroconvert
  vector<lower = 0, upper = 1>[n_ab > 0 ? 2 : 0] gamma; // antibody waning (inf & vac)
  vector<lower = 0, upper = 1>[n_ab > 0 ? 1 : 0] delta; // vaccine efficacy
  vector<lower = 0, upper = 1>[n] init_dab; // initial proportion with antibodies
  vector<lower = 0>[n_ab > 0 ? 1 : 0] k; // Potential loss of efficacy from new infections in already seropositive people
  vector<lower = 0>[n_ab > 0 ? 1 : 0] l; // Potential loss of efficacy from new doses being administered to already seropositive people
}

transformed parameters {
  matrix[n, t] gp; // value of gp at time t + initialisation
  matrix[n, t] infections; // incident infections at time t
  matrix[n_ab, t] infs_with_potential_abs; // Infections with the potential to have ab
  matrix[n, t] dcases; // detectable cases at time t
  matrix[n_ab, t] dab; // proportion of individuals with antibodies at time t
  matrix[n, obs] odcases;
  matrix[n_ab, ab_obs] odab;
  matrix[n, obs] combined_sigma;
  matrix[n_ab, ab_obs] combined_ab_sigma;
  // update gaussian process
  for (i in 1:n) {
    gp[i, (1 + diff_order):t] = to_row_vector(update_gp(PHI, M, L, alpha[i], rho[i], to_vector(eta[i]), 0));
    // setup differencing of the GP
    if (diff_order) {
      gp[i, 1:diff_order] = to_row_vector(init_growth[i]);
      for (j in 1:diff_order) {
        gp[i] = cumulative_sum(gp[i]);
      }
    }

    // relative probability of infection
    // inc_init is the mean incidence
    infections[i] = inv_logit(init_inc[i] + gp[i]);

    // calculate detectable cases
    dcases[i] = to_row_vector(convolve(to_vector(infections[i]), prob_detect));
    // calculate observed detectable cases
    odcases[i] = to_row_vector(observed_in_window(to_vector(dcases[i]), prev_stime, prev_etime, ut, obs));
    //combined standard error
    combined_sigma[i] = sqrt(rep_row_vector(square(sigma), obs) + prev_sd2[i]);
  }

  //calculate infections with potential to have antibodies
  for (i in 1:n_ab) {
    infs_with_potential_abs[i] = to_row_vector(convolve(to_vector(infections[ab_index[i]]), inf_ab_delay));
    // calculate detectable antibodies
    dab[i] = to_row_vector(detectable_antibodies(to_vector(infs_with_potential_abs[i]),
                                                 to_vector(vacc_with_ab[ab_index[i]]), 
                                                 beta[1], gamma, delta[1], k[1], l[1],
                                                 init_dab[ab_index[i]], t));
    // calculate observed detectable antibodies
    odab[i] = to_row_vector(observed_in_window(to_vector(dab[i]), ab_stime, ab_etime, ut, ab_obs));
    //combined standard error
    combined_ab_sigma[i] = sqrt(rep_row_vector(square(ab_sigma[1]), ab_obs) + ab_sd2[i]);
  }
}

model {
  // gaussian process priors
  for (i in 1:n) {
    rho[i] ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
    alpha[i] ~ std_normal() T[0,];
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
    logit(beta) ~ normal(pbeta[1], pbeta[2]);
    logit(gamma) ~ normal(pgamma_mean, pgamma_sd); 
    logit(delta) ~ normal(pdelta[1], pdelta[2]);
    log(k) ~ normal(0, 0.1);
    log(l) ~ normal(0, 0.1);
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
  matrix[n, t - ut] R;
  matrix[n, t - ut] r;
  matrix[n, t] cumulative_infections;
  real est_prev[n, obs];
  real est_ab[n, ab_obs];

  {
    // sample generation time
    real gtm_sample = normal_rng(gtm[1], gtm[2]);
    real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);

    // get cumulative incidence
    for (i in 1:n) {
      cumulative_infections[i] = cumulative_sum(infections[i]);
      // sample estimated prevalence
      est_prev[i] = normal_rng(odcases[i], sigma);
      // calculate Rt using infections and generation time
      R[i] = to_row_vector(calculate_Rt(to_vector(infections[i]), ut, gtm_sample, gtsd_sample, gtmax, 1));
      // calculate growth
      r[i] = to_row_vector(calculate_growth(to_vector(infections[i]), ut));
    }

    if (n_ab > 0) {
      matrix[n, t] gen_infs_with_potential_abs; // Infections with the potential to have ab
      matrix[n, t] gen_dab;
      matrix[n, ab_obs] gen_odab;

      for (i in 1:n) {
        gen_infs_with_potential_abs[i] = to_row_vector(convolve(to_vector(infections[i]), inf_ab_delay));
        // calculate detectable antibodies
        gen_dab[i] = to_row_vector(detectable_antibodies(to_vector(gen_infs_with_potential_abs[i]), 
                                                         to_vector(vacc_with_ab[i]), beta[1],
                                                         gamma, delta[1], k[1], l[1], init_dab[i], t));
        // calculate observed detectable antibodies
        gen_odab[i] = to_row_vector(observed_in_window(to_vector(gen_dab[i]), ab_stime, ab_etime, ut, ab_obs));
        // estimate with error
        est_ab[i] = normal_rng(gen_odab[i], ab_sigma[1]);
      }
    }
  }
}
