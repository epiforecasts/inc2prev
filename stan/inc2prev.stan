functions {
#include functions/gaussian_process.stan
#include functions/rt.stan
#include functions/prev.stan
#include functions/generated_quantities.stan
}

data {
#include data/observations.stan
  vector[obs] prev;
#include data/observation_generation.stan
#include data/summary_measures.stan
#includedata/gaussian_process.stan
}
  
transformed data {
#include transformed-data/gaussian_process.stan
}

parameters {
#include parameters/gaussian_process.stan
#include parameters/observation_model.stan
}

transformed parameters {
#include transformed-parameters/guassian_process.stan
#include transformed-parameters/observation_generation.stan
}

model {
#include model/gaussian_process.stan
#include model/observation_model.stan
}

generated quantities {
#include generated-quantities/summary_measures.stan
}
