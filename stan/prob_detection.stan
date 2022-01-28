functions {
#include functions/prob_detection.stan
}
data {
#include chunks/data/prob_detection.stan
}

parameters {
  real inc_mean;
  real<lower = 0> inc_sd;
  vector[3] effs;
  real<lower = 0> change;
  vector <lower = 0, upper = 1> [pb_p] inf_at; 
}

model {
  detection_prob_lp(
    // Observations
    pb_result, pb_test_day, pb_sym_at_test,
    pb_last_asym_at_test, pb_inf_upper_bound, 
    // Priors
    inc_mean, inc_sd, inf_at, effs, change,
    // Indexs
    pb_p, pb_n, pb_id,
    // Prior parameterisation
    inc_mean_p, inc_sd_p
  );
}

generated quantities {
  vector[301] pb;
  pb = detection_prob(30, 10, 0, effs, change);
}
