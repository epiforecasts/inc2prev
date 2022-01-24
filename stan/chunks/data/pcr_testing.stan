  int pcr_p; // Number of study members
  int pcr_n; // data points in study
  int pcr_id[pcr_n]; //study ids
  vector[pcr_n] pcr_test_day; // day of test
  int pcr_result[pcr_n]; // Result of PCR test
  vector[pcr_p] pcr_sym_at_test; // day of first symptomatic test
  vector[pcr_p] pcr_last_asym_at_test; // day of last test with no symptoms (currently or previously)
  vector[pcr_p] pcr_inf_upper_bound;  // maximum time at which infection must have occured before
  vector[2] inc_mean_p; // Mean and sd of the incubation period mean
  vector[2] inc_sd_p; // Mean and sd of the incubation period mean