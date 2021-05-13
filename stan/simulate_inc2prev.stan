functions {
#include functions/rt.stan
#include functions/prev.stan
#include functions/generated_quantities.stan
}

data {
#include data/observations.sta
#include data/observation_generation.stan
#include data/summary_measures.stan
#include data/gaussian_process.stan
  vector[t] gp;
}
  
generated quantities {
#include transformed-parameters/observation_generation.stan 
#include generated-quantities/summary_measures.stan
#include generated-quantities/observation.stan
}
