functions {
#include gaussian_process.stan
#include rt.stan
#include generated_quantities.stan
}

data {
  int ut;
  int ot;
  int t;
  int obs;
  real prev[obs];
  int prev_time[obs];
  int pbt;
  real prob_detect[pbt + 1];
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
  }
  
transformed data {
  // set up approximate gaussian process
  matrix[t, M] PHI = setup_gp(M, L, t);  
}
  

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
}

transformed parameters {
  vector[t] gp;
  vector[t] prob_inf;
  // update gaussian process
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // relative probability of infection
  prob_inf = inv_logit(gp);
  prob_inf = prob_inf ./ sum(prob_inf);
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
}

generated quantities {
  vector[t - 7] R;
  vector[t - 1] r;
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(prob_inf, 7, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(prob_inf, 1);
}
