functions {
#include functions/detection_prob.stan
}
data {
#include chunks/data/pcr_testing.stan
}

parameters {
  real inc_mean;
  real<lower = 0> inc_sd;
  vector[3] effs;
  real<lower = 0> change;
  vector <lower = 0, upper = 1> [pcr_p] inf_at; 
}

model {
  detection_prob_lp(
    // Observations
    pcr_result, pcr_test_day, pcr_sym_at_test,
    pcr_last_asym_at_test, pcr_inf_upper_bound, 
    // Priors
    inc_mean, inc_sd, inf_at, effs, change,
    // Indexs
    pcr_p, pcr_n, pcr_id,
    // Prior parameterisation
    inc_mean_p, inc_sd_p
  );
}

generated quantities {
  vector[301] pb;
  pb = detection_prob(30, 10, 0, effs, change);
}
