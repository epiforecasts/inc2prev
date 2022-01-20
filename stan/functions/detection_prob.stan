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

vector detection_prob_by_day(int days, vector effs, real bp) {
  vector[days+1] pb;
  vector[days+1] k;
  for (i in 0:(days)) {
    k[i+1] = i + 0.5; // Probability at halfway point of the day
  }
  pb = detection_prob_logit(k, effs, bp);
  pb = inv_logit(pb);
  return(pb);
}

