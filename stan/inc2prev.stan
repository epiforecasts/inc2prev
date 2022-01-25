functions {
#include detection_prob.stan
#include rev_vec.stan
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
  int obs; // number of prevalence observations
  int ab_obs; // number of antibody prevalence observations
  vector[obs] prev; // observed positivity prevalence
  vector[obs] prev_sd2; // squared standard deviation of observed positivity prevalence
  vector[ab_obs] ab; // observed antibody posivitiy prevalence
  vector[ab_obs] ab_sd2; // squared standard deviation of observed antibody prevalence
  int prev_stime[obs]; // starting times of positivity prevalence observations
  int prev_etime[obs]; // end times of positivity prevalence observations
  int ab_stime[ab_obs]; // starting times of antibody prevalence observations
  int ab_etime[ab_obs]; // end times of antibody prevalence observations
  vector[ab_obs] vacc; // vaccinations
  int pbt; // maximum detection time
  vector[3] pcr_eff_m; //Mean detection probability effects
  vector[3] pcr_eff_sd; //SD detection probability effects
  real pcr_change_m; // breakpoint of PCR detection
  real pcr_change_sd; // breakpoint of PCR detection
  vector[pbt+1] prob_detect_mean; // at each time since infection, probability of detection
  vector[pbt+1] prob_detect_sd; // at each time since infection, tandard deviation of probability of detection
  // Mode to use for probability of detection, 0 = exact,
  // 1 = posterior summary, and 2 = from study data
  int pb_mode; 
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  int diff_order; // Order of differencing to use for the GP (0 = cases -> mean, 1 = growth_t -> 0, 2 = growth_t -> growth_{t-1})
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
  real init_inc_mean ; // Mean initial/mean incidence (logit)
  real init_ab_mean[ab_obs ? 1 : 0]; // mean estimate of initial antibody prevalence
  real init_ab_sd[ab_obs ? 1 : 0]; // sd of estimate of initial antibody prevalence
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
  vector[ab_obs ? t : 0] vacc_with_ab;
  // set up approximate gaussian process
  matrix[t - diff_order, M] PHI = setup_gp(M, L, t - diff_order);
  // Calculate vaccinations with the potential to have antibodies
  if (ab_obs) {
    vacc_with_ab = convolve(vacc, vacc_ab_delay);
  }
}

parameters {
  real<lower = 0> rho; // length scale of gp
  real<lower = 0> alpha; // scale of gp
  vector[M] eta; // eta of gp
  real init_inc; // Initial infections
  vector[diff_order] init_growth;
  real<lower = 0> sigma; // observation error
  real<lower = 0> pb_sigma;
  vector<lower = 0>[ab_obs ? 1 : 0] ab_sigma; // observation error
  vector[3] pcr_eff; // probability of detection piecewise linear effects
  real<lower = 0> pcr_change; // probability of detection breakpoint
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] beta; // proportion that don't seroconvert
  vector<lower = 0, upper = 1>[ab_obs ? 2 : 0] gamma; // antibody waning (inf & vac)
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] delta; // vaccine efficacy
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] init_dab; // initial proportion with antibodies
}

transformed parameters {
  vector[t] gp; // value of gp at time t + initialisation 
  vector[t] infections; // incident infections at time t
  vector[ab_obs ? t : 0] infs_with_potential_abs; // Infections with the potential to have ab
  vector<lower = 0, upper = 1>[pbt+1] prob_detect;
  vector[pbt+1] combined_pb_sigma;
  vector<lower = 0, upper = 1>[t] dcases; // detectable cases at time t
  vector[ab_obs ? t : 0] dab; // proportion of individuals with antibodies at time t
  vector[obs] odcases;
  vector[ab_obs] odab;
  vector[obs] combined_sigma;
  vector[ab_obs] combined_ab_sigma;
  // update gaussian process
  gp[(1 + diff_order):t] = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // setup differencing of the GP
  if (diff_order) {
    gp[1:diff_order] = init_growth;
    for (i in 1:diff_order) {
      gp = cumulative_sum(gp);
    }
  }

  // relative probability of infection
  // inc_init is the mean incidence
  infections = inv_logit(init_inc + gp);
  // calculate probability of detection
  prob_detect = detection_prob_by_day(pbt, pcr_eff, pcr_change);
  combined_pb_sigma = sqrt(square(pb_sigma) + square(prob_detect_sd));
  // calculate detectable cases
  dcases = convolve(infections, rev_vec(prob_detect));
  // calculate observed detectable cases
  odcases = observed_in_window(dcases, prev_stime, prev_etime, ut, obs);
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);

  //calculate infections with potential to have antibodies
  if (ab_obs) {
    infs_with_potential_abs = convolve(infections, inf_ab_delay);
    // calculate detectable antibodies
    dab = detectable_antibodies(infs_with_potential_abs, vacc_with_ab, beta[1],
                                gamma, delta[1], init_dab[1], t);
    // calculate observed detectable antibodies
    odab = observed_in_window(dab, ab_stime, ab_etime, ut, ab_obs);
    //combined standard error
    combined_ab_sigma = sqrt(square(ab_sigma[1]) + ab_sd2);
  }
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();

  // Initial infections
  init_inc ~ normal(init_inc_mean, 2);
  if (diff_order) {
    init_growth ~ normal(0, 0.25);
  }

  // prevalence observation model
  pcr_eff ~ normal(pcr_eff_m, pcr_eff_sd);
  pcr_change ~ normal(pcr_change_m, pcr_change_sd);
  pb_sigma ~ normal(0.025, 0.025) T[0, ];
  prob_detect_mean ~ normal(prob_detect, combined_pb_sigma);

  // Priors for antibody model
  if (ab_obs) {
    init_dab ~ normal(init_ab_mean[1], init_ab_sd[1]);
    logit(beta) ~ normal(pbeta[1], pbeta[2]);
    logit(gamma) ~ normal(pgamma_mean, pgamma_sd); 
    logit(delta) ~ normal(pdelta, pdelta); 
    ab_sigma[1] ~ normal(0.025, 0.025) T[0,];
  }

  sigma ~ normal(0.005, 0.0025) T[0,];
 
  if (prev_likelihood) {
    prev ~ normal(odcases, combined_sigma);
  }
  if (ab_likelihood && ab_obs) {
    ab ~ normal(odab, combined_ab_sigma);
  }
}

generated quantities {
  vector[t - ut] R;
  vector[t - ut] r;
  vector[t] cumulative_infections;
  real est_prev[obs];
  real est_ab[ab_obs];
  // get cumulative incidence
  cumulative_infections = cumulative_sum(infections);
  // sample estimated prevalence
  est_prev = normal_rng(odcases, combined_sigma);
  if (ab_obs) {
    est_ab = normal_rng(odab, combined_ab_sigma);
  }
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(infections, ut, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(infections, ut);
}
