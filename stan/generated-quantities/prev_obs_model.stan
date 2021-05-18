  // prevalence observation model
  real sigma;
  vector[obs] prev;
#include tparameters-var-def/prev_obs_model.stan 
  sigma = normal_rng(0.005, 0.0025) T[0,];
#include tparameters/prev_obs_model.stan 
  prev = normal_rng(odcases, combined_sigma);
