  // prevalence observation model
  sigma ~ normal(0.005, 0.0025) T[0,];
  prev ~ normal(odcases, combined_sigma) T[0, 1];