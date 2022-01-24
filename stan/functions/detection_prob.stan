// Probability of detection using a piecewise linear with a single breakpoint
// on a logit scale
// Based on: https://doi.org/10.1186/s12916-021-01982-x
vector detection_prob_logit(vector k, vector effs, real bp) {
  int l = num_elements(k);
  vector[l] t;
  vector[l] pb;
  t = k - bp; // centre on breakpoint
  for (i in 1:l) {
    pb[i] = effs[1] + effs[2] * t[i] + t[i] * effs[3] * effs[2] * step(t[i]);
  }
  return(pb);
}

vector detection_prob(int l, int interval, real shift, vector effs, real bp) {
  int ind = l*interval;
  vector[ind+1] pb;
  vector[ind+1] k;
  for (i in 0:ind) {
    k[i+1] = (i + shift);
  }
  k = k / interval;
  pb = detection_prob_logit(k, effs, bp);
  pb = inv_logit(pb);
  return(pb);
}


vector detection_prob_by_day(int days, vector effs, real bp) {
  vector[days+1] pb;
  // detection probability at the half point of each day
  pb = detection_prob(days, 1, 0.5, effs, bp);
  return(pb);
}

void detection_prob_lp(int[] result, vector test_day, vector sym_at_test,
                       vector last_asym_at_test, vector inf_upper_bound, 
                       real inc_mean, real inc_sd, vector inf_at, 
                       vector effs, real change, int p, int n, int[] id,
                       vector inc_mean_p, vector inc_sd_p) {

  // Priors on the incubation period
  inc_mean ~ normal(inc_mean_p[1], inc_mean_p[2]);
  inc_sd ~ normal(inc_sd_p[1], inc_sd_p[2]) T[0, ];

  // Prior on time infected at (as a proportion of time
  // from start to upper bound)
  inf_at ~ beta(3, 1);

  // Priors on piecewise linear (on logit) probability
  effs ~ std_normal();
  change ~ normal(5, 5) T[0, ];

  // Infection time
  vector[p] tinf = inf_upper_bound .* inf_at;

  // Symptom onset likelihood
  // Probability of onset before first symptoms minus probability prior to
  // last asymptomatic test (or 0 if occurred prior to infection).
  vector[p] inf_to_sym = sym_at_test - tinf;
  vector[p] inf_to_lasym = last_asym_at_test - tinf;
  for (i in 1:p) {
    target += log(
      lognormal_cdf(inf_to_sym[i] | inc_mean, inc_sd) - 
      (inf_to_lasym[i] <= 0 ? 
        0 : lognormal_cdf(inf_to_lasym[i] | inc_mean,  inc_sd)
      )
    );
  }

  // detection likelihood
  vector[n] inf_to_test = test_day - tinf[id];
  vector[n] p_d = detection_prob_logit(inf_to_test, effs, change);
  result ~ bernoulli_logit(p_d);
  // Add negative test at date of infection to constrain PCR detection
  real test_at_zero = effs[1] - change * effs[2];
  target += n * bernoulli_logit_lpmf(0 | test_at_zero);
}