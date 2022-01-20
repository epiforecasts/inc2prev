  int pcr_n; // Number of study members
  int pcr_p; // data points in study
  int pcr_id[pcr_p ? pcr_p : 1]; //study ids
  vector[pcr_p] pcr_test_day; // day of test
  int pcr_result[pcr_p]; // Result of PCR test
  vector[pcr_n ? pcr_n : 1] pcr_sym_at_test; // day of first symptomatic test
  vector[pcr_n ? pcr_n : 1] pcr_last_asym_at_test; // day of last test with no symptoms (currently or previously)
  vector[pcr_n ? pcr_n : 1] pcr_inf_upper_bound;  // maximum time at which infection must have occured before
  real inc_mean_p[2]; // Mean and sd of the incubation period mean
  real inc_sd_p[2]; // Mean and sd of the incubation period mean