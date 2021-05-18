  // relative probability of infection
  infections = N * inv_logit(inc_zero + covariates);
  // calculate detectable cases
  dcases = detectable_cases(infections, prob_detect, pbt, t);
  // calculate observed detectable cases
  odcases = observed_cases(dcases, prev_stime, prev_etime, ut, obs);
  odcases = odcases / N;
