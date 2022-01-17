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
  int c_obs; // Number of count observations
  int c_grps; // Number of groups of counts
  int c_obs_stime[c_grps]; // Group starting membership for each count
  int c_obs_etime[c_grps]; // Group ending membership for each count
  int c_time[c_obs]; // Occurance times for counts
  vector[c_obs] counts;
  vector[t] vacc; // vaccinations
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
  real init_inc_mean ; // Mean initial/mean incidence (logit)
  real init_ab_mean; // mean estimate of initial antibody prevalence
  real init_ab_sd;   // sd of estimate of initial antibody prevalence
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
  int count_likelihood; // Should the likelihood for count data be included
}

transformed data {
  vector[t] vacc_with_ab;
  // set up approximate gaussian process
  matrix[t - diff_order, M] PHI = setup_gp(M, L, t - diff_order);
  // Calculate vaccinations with the potential to have antibodies
  vacc_with_ab = convolve(vacc, vacc_ab_delay);
}

parameters {
  real<lower = 0> rho; // length scale of gp
  real<lower = 0> alpha; // scale of gp
  vector[M] eta; // eta of gp
  real init_inc; // Initial infections
  vector[diff_order] init_growth;
  real<lower = 0> sigma; // observation error
  real<lower = 0> ab_sigma; // observation error
  vector<lower = 0, upper = 1>[pbt] prob_detect; // probability of detection as a function of time since infection
  real<lower = 0, upper = 1> beta; // proportion that don't seroconvert
  vector<lower = 0, upper = 1>[2] gamma; // antibody waning (inf & vac)
  real<lower = 0, upper = 1> delta; // vaccine efficacy
  real<lower = 0, upper = 1> init_dab; // initial proportion with antibodies
}

transformed parameters {
  vector[t] inf_gp; // value of gp at time t + initialisation 
  vector[t] cscale[c_grps]; //time-varying scaling factor of infections to counts
  vector[t] infections; // incident infections at time t
  vector[t] infs_with_potential_abs; // Infections with the potential to have ab
  vector<lower = 0, upper = 1>[t] dcases; // detectable cases at time t
  vector[t] dab; // proportion of individuals with antibodies at time t
  vector[obs] odcases;
  vector[ab_obs] odab;
  vector[obs] combined_sigma;
  vector[ab_obs] combined_ab_sigma;
  // update gaussian process
  inf_gp[(1 + diff_order):t] = update_gp(PHI, M, L, alpha[1], rho[1], eta[1:M],
                                         0);
  // setup differencing of the GP
  if (diff_order) {
    inf_gp[1:diff_order] = init_growth;
    for (i in 1:diff_order) {
      inf_gp = cumulative_sum(inf_gp);
    }
  }
  // relative probability of infection
  // inc_init is the mean incidence
  infections = inv_logit(init_inc + inf_gp);
  // calculate detectable cases
  dcases = convolve(infections, prob_detect);
  // calculate observed detectable cases
  odcases = observed_in_window(dcases, prev_stime, prev_etime, ut, obs);
  //calculate infections with potential to have antibodies
  infs_with_potential_abs = convolve(infections, inf_ab_delay);
  // calculate detectable antibodies
  dab = detectable_antibodies(infs_with_potential_abs, vacc_with_ab, beta,
                              gamma, delta, init_dab, t);
  // calculate observed detectable antibodies
  odab = observed_in_window(dab, ab_stime, ab_etime, ut, ab_obs);
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);
  combined_ab_sigma = sqrt(square(ab_sigma) + ab_sd2);
  // calculate infections reported as counts
  if (c_grps) {
    for (i in 1:c_grps) {
      cscale[i] = update_gp(cPHI[c_grps], M, L, alpha[1+c_grps],
                            rho[1+c_grps],
                            eta[((c_grps-1)*M + 1):(c_grps*M)], 0);
      cscale[i] = inv_logit(cscale_mean[i] + cscale[i]);
      cdist = discretised_gamma_pmf(cdist_indexes, cdist_mean[i], cdist_sd[i],
                                    cdist_max, 0);
      ecounts[c_grps] = convolve(cscale[c_grps] .* infections, cdist[c_grps]);
    }
  }
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  for (i in 1:gps) {
    alpha[i] ~ std_normal() T[0,];
  }
  eta ~ std_normal();

  // Initial infections
  init_inc ~ normal(init_inc_mean, 2);
  if (diff_order) {
    init_growth ~ normal(0, 0.25);
  }
  // prevalence observation model
  for (i in 1:pbt) {
    prob_detect[i] ~ normal(prob_detect_mean[i], prob_detect_sd[i]) T[0, 1];
  }

  // Priors for antibody model
  init_dab ~ normal(init_ab_mean, init_ab_sd);
  logit(beta) ~ normal(pbeta[1], pbeta[2]);
  logit(gamma) ~ normal(pgamma_mean, pgamma_sd); 
  logit(delta) ~ normal(pdelta, pdelta); 


  // Priors for count data model
  if (c_grps) {
    cscale_mean ~ normal(0, 5);
    for (i in c_grps) {
      cdist_mean[i] ~ normal(5, 10) T[0,];
      cdist_sd[i] ~ normal(5, 10) T[0,];
    }
  }

  sigma ~ normal(0.005, 0.0025) T[0,];
  ab_sigma ~ normal(0.025, 0.025) T[0,];
  if (prev_likelihood) {
    prev ~ normal(odcases, combined_sigma);
  }
  if (ab_likelihood) {
    ab ~ normal(odab, combined_ab_sigma);
  }
  if (count_likelihood & c_grps) {
    for (i in 1:c_grps) {
      counts[c_obs_stime[i]:c_obs_etime[i]] ~ neg_binomial_2(ecounts[i], phi[i]);
    }
  }
}

generated quantities {
  vector[t - ut] R;
  vector[t - ut] r;
  real est_prev[obs];
  real est_ab[ab_obs];
  // sample estimated prevalence
  est_prev = normal_rng(odcases, combined_sigma);
  est_ab = normal_rng(odab, combined_ab_sigma);
  if (c_grps) {
    for (i in 1:c_grps) {
      est_counts[c_obs_stime[i]:c_obs_etime[i]] ~ neg_binomial_2(ecounts[i], phi[i]);
        }
  }
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(infections, ut, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(infections, ut);
}
