functions {
#include functions/gaussian_process.stan
#include functions/rt.stan
#include functions/prev.stan
#include functions/generated_quantities.stan
}

data {
#include data/observations.stan
#include data/prev_obs_model.stan
#include data/observation_generation.stan
#include data/summary_measures.stan
#include data/gaussian_process.stan
}
  
transformed data {
#include transformed-data/gaussian_process.stan
}

parameters {
#include parameters/gaussian_process.stan
#include parameters/prev_obs_model.stan
}

transformed parameters {
#include tparameters-var-def/gaussian_process.stan
#include tparameters-var-def/observation_generation.stan
#include tparameters-var-def/prev_obs_model.stan
#include tparameters/gaussian_process.stan
#include tparameters/observation_generation.stan
#include tparameters/prev_obs_model.stan
}

model {
#include model/gaussian_process.stan
#include model/prev_obs_model.stan
}

generated quantities {
#include generated-quantities-var-def/summary_measures.stan
#include generated-quantities/summary_measures.stan
}
