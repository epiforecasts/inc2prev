  vector[t] infections;
  vector[t] dcases;
  vector[obs] odcases;
  vector[obs] combined_sigma;
  // relative probability of infection
  infections = N * inv_logit(inc_zero + gp);
  // calculate detectable cases
  dcases = detectable_cases(infections, prob_detect, pbt, t);
  // calculate observed detectable cases
  odcases = observed_cases(dcases, prev_stime, prev_etime, ut, obs);
  odcases = odcases / N;
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);
