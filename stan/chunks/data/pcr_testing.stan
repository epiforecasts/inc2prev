  int pb_p; // Number of study members
  int pb_n; // data points in study
  int pb_id[pb_n]; //study ids
  vector[pb_n] pb_test_day; // day of test
  int pb_result[pb_n]; // Result of pb test
  vector[pb_p] pb_sym_at_test; // day of first symptomatic test
  vector[pb_p] pb_last_asym_at_test; // day of last test with no symptoms (currently or previously)
  vector[pb_p] pb_inf_upper_bound;  // maximum time at which infection must have occured before
  vector[2] inc_mean_p; // Mean and sd of the incubation period mean
  vector[2] inc_sd_p; // Mean and sd of the incubation period mean