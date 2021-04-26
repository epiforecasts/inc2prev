functions {
#include gaussian_process.stan
#include rt.stan
#include prev.stan
#include generated_quantities.stan
}

data {
  int ut;
  int ot;
  int t;
  int obs;
  vector[obs] prev;
  vector[obs] prev_sd2;
  int prev_stime[obs];
  int prev_etime[obs];
  int pbt;
  vector[pbt] prob_detect;
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
  int N;
  real inc_zero;
}
  
transformed data {
  // set up approximate gaussian process
  matrix[t, M] PHI = setup_gp(M, L, t); 
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  real<lower = 0> sigma;
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
  infections = N * inv_logit(inc_zero + gp);
  // calculate detectable cases
  dcases = detectable_cases(infections, prob_detect, pbt, t);
  // calculate observed detectable cases
  odcases = observed_cases(dcases, prev_stime, prev_etime, ut, obs);
  odcases = odcases / N;
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();

  // prevalence observation model
  sigma ~ normal(0.005, 0.0025) T[0,];
  prev ~ normal(odcases, combined_sigma);
}

generated quantities {
  vector[t - 7] R;
  vector[t - 1] r;
  vector[obs] est_prev;
  vector[ot] pop_prev;
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(infections, 7, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(infections, 1);
  // population prevelence
  pop_prev = dcases[(ut + 1):t] / N;
  // sample estimated prevalence
  est_prev = to_vector(normal_rng(odcases, combined_sigma));
}
