  // sample generation time
  gtm_sample = normal_rng(gtm[1], gtm[2]);
  gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(infections, 7, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(infections, 1);
  // population prevelence
  pop_prev = dcases[(ut + 1):t] / N;
  // sample estimated prevalence
  est_prev = to_vector(normal_rng(odcases, combined_sigma));
  // population incidence rate (per 100,000k)
  pop_inc = infections[(ut + 1):t] / N * 100000;
