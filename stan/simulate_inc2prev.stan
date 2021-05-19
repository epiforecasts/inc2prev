functions {
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
#include tparameters-var-def/gaussian_process.stan
}
  
generated quantities {
#include tparameters-var-def/observation_generation.stan 
#include generated-quantities-var-def/summary_measures.stan
#include generated-quantities-var-def/prev_obs_model.stan
#include parameters/prev_obs_model.stan
#include tparameters/observation_generation.stan  
#include generated-quantities/prev_obs_model.stan 
#include generated-quantities/summary_measures.stan
}
