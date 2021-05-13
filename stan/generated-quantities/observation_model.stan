  // prevalence observation model
  real sigma;
  vector[obs] prev;
  sigma = normal_rng(0.005, 0.0025) T[0,];
  prev = normal_rng(odcases, combined_sigma);
