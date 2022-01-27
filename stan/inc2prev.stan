functions {
#include prob_detection.stan
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
  vector[3] pb_effs_m; //Mean detection probability effects
  vector[3] pb_effs_sd; //SD detection probability effects
  real pb_change_m; // breakpoint of pb detection
  real pb_change_sd; // breakpoint of pb detection
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
  // Latent expected infections
  real<lower = 0> rho; // length scale of gp
  real<lower = 0> alpha; // scale of gp
  vector[M] eta; // eta of gp
  real init_inc; // Initial infections
  vector[diff_order] init_growth;

  
  // cross sectional prevalence
  vector[3] pb_effs; // probability of detection piecewise linear effects
  real<lower = 0> pb_change; // probability of detection breakpoint
  real<lower = 0> pb_sigma;
  real<lower = 0> prev_sigma; // observation error
  vector[pb_mode == 2 ? 1 : 0] inc_mean; // Incubation period mean
  real<lower = 0>[pb_mode == 2 ? 1 : 0] inc_sd; // Incubation period SD
  vector <lower = 0, upper = 1> [pb_p] inf_at; 

  // Antibody model parameters
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] beta; // proportion that don't seroconvert
  vector<lower = 0, upper = 1>[ab_obs ? 2 : 0] gamma; // antibody waning (inf & vac)
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] delta; // vaccine efficacy
  vector<lower = 0, upper = 1>[ab_obs ? 1 : 0] init_dab; // initial proportion with antibodies
  vector<lower = 0>[ab_obs ? 1 : 0] ab_sigma; // observation error
}

transformed parameters {
  vector[t] gp; // value of gp at time t + initialisation 
  vector[t] infections; // incident infections at time t
  vector[ab_obs ? t : 0] infs_with_potential_abs; // Infections with the potential to have ab
  vector<lower = 0, upper = 1>[pbt+1] prob_detect;
  vector[pb_mode == 1 ? pbt+1 : 0] combined_pb_sigma;
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
  if (pb_mode) {
    prob_detect = prob_detect_mean;
  }else{
    prob_detect = detection_prob_by_day(pbt, pb_effs, pb_change);
  }
  
  if (pb_mode == 1) {
    combined_pb_sigma = sqrt(square(pb_sigma) + square(prob_detect_sd));
  }
  
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
    //combined standard errors
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

  // prevalence detection prob model
  if (pb_mode == 1) {
    // summary data
    pb_effs ~ normal(pb_effs_m, pb_effs_sd);
    pb_change ~ normal(pb_change_m, pb_change_sd);
    pb_sigma ~ normal(0.025, 0.025) T[0, ];
    prob_detect_mean ~ normal(prob_detect, combined_pb_sigma);
  }else if (pb_mode == 2) {
    // individual data
    detection_prob_lp(
      // Observations
      pb_result, pb_test_day, pb_sym_at_test,
      pb_last_asym_at_test, pb_inf_upper_bound, 
      // Priors
      inc_mean, inc_sd, inf_at, pb_effs, pb_change,
      // Indexs
      pb_p, pb_n, pb_id,
      // Prior parameterisation
      inc_mean_p, inc_sd_p
    );
  }    
  // prevalence observation model
  prev_sigma ~ normal(0.005, 0.0025) T[0,];
 

  // Priors for antibody model
  if (ab_obs) {
    init_dab ~ normal(init_ab_mean[1], init_ab_sd[1]);
    logit(beta) ~ normal(pbeta[1], pbeta[2]);
    logit(gamma) ~ normal(pgamma_mean, pgamma_sd); 
    logit(delta) ~ normal(pdelta, pdelta); 
    ab_sigma[1] ~ normal(0.025, 0.025) T[0,];
  }

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
