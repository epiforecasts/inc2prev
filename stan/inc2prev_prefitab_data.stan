functions {
#include gaussian_process.stan
#include rt.stan
#include convolve.stan
#include observed_in_window.stan
#include ab.stan
#include generated_quantities.stan
}

data {
  int ut;
  int t;
  int obs;
  int ab_obs; // number of antibody prevalence observations
  vector[obs] prev;
  vector[obs] prev_sd2;
  int prev_stime[obs];
  int prev_etime[obs];
  int ab_stime[ab_obs]; // starting times of antibody prevalence observations
  int ab_etime[ab_obs]; // end times of antibody prevalence observations
  int pbt;
  vector[pbt] prob_detect_mean;
  vector[pbt] prob_detect_sd;
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
  real inc_zero;

    
  int n_ab_draws;
  int n_ab_pars;
  matrix[n_ab_draws, n_ab_pars] ab_draws;

  
  
  vector[t] vacc; // vaccinations
  int linf_ab_delay; // Length of the delay between infection and ab
  // PMF of the delay between infection and ab
  vector[linf_ab_delay] inf_ab_delay; 
  int lvacc_ab_delay; // Length of the delay between vaccination and ab
  // PMF of the  delay between vaccination and ab
  vector[lvacc_ab_delay] vacc_ab_delay; 
  
}

transformed data {
  vector[t] vacc_with_ab;
  // set up approximate gaussian process
  matrix[t, M] PHI = setup_gp(M, L, t);
    // Calculate vaccinations with the potential to have antibodies
  vacc_with_ab = convolve(vacc, vacc_ab_delay);

}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  real<lower = 0> sigma;
  vector<lower = 0, upper = 1>[pbt] prob_detect;
}

transformed parameters {
  vector[t] gp;
  vector[t] infections;
  vector[t] dcases;
  vector[obs] odcases;
  vector[obs] combined_sigma;
  // update gaussian process
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // relative probability of infection
  infections = inv_logit(inc_zero + gp);
  // calculate detectable cases
  dcases = convolve(infections, prob_detect);
  // calculate observed detectable cases
  odcases = observed_in_window(dcases, prev_stime, prev_etime, ut, obs);
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);
  
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();

  // prevalence observation model
  for (i in 1:pbt) {
    prob_detect[i] ~ normal(prob_detect_mean[i], prob_detect_sd[i]) T[0, 1];
  }
  sigma ~ normal(0.005, 0.0025) T[0,];
  prev ~ normal(odcases, combined_sigma);
}

generated quantities {
  vector[t - ut] R;
  vector[t - 1] r;
  vector<lower = 0>[t]  cumulative_infections;
  real est_prev[obs];
  
  vector[t] infs_with_potential_abs; // Infections with the potential to have ab
  vector[t] dab; // proportion of individuals with antibodies at time t
  vector[ab_obs] odab;
  vector[ab_obs] combined_ab_sigma;
  real init_dab;
  
  real ab_index_r;
  int ab_index;
  real beta; 
  vector[2] gamma; 
  real delta;
  real ab_sigma; 
  
  real est_ab[ab_obs];
  
  
  // cumulative incidence
  cumulative_infections = cumulative_sum(infections);

  // sample estimated prevalence
  est_prev = normal_rng(odcases, combined_sigma);
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(infections, ut, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(infections, 1);
  
  ab_index_r = n_ab_draws * uniform_rng(0,1);
  ab_index = 0;
  while(ab_index < ab_index_r)
    ab_index += 1;
  
  ab_index += 1;
  
  beta = ab_draws[ab_index][1];
  gamma = to_vector(ab_draws[ab_index][2:3]);
  delta = ab_draws[ab_index][4];
  ab_sigma = ab_draws[ab_index][5];
  init_dab = ab_draws[ab_index][6];
    
  infs_with_potential_abs = convolve(infections, inf_ab_delay);
  // calculate detectable antibodies
  dab = detectable_antibodies(infs_with_potential_abs, vacc_with_ab, beta,
                             gamma, delta, init_dab, t);
  // calculate observed detectable antibodies
  odab = observed_in_window(dab, ab_stime, ab_etime, ut, ab_obs);
  //combined standard error
  //combined_ab_sigma = ab_sigma;
  
  est_ab = normal_rng(odab, ab_sigma);
  
}
