// Antibodies based on cumulative infections and vaccination
vector detectable_antibodies(vector infections, real[] vacc,
                             real beta, real gamma, real delta, real init_pop_ab, int t) {
  vector[t] pop_ab;

  pop_ab[1] =
    init_pop_ab // initial antibodies
    + (1 - beta) * infections[1] / (1 - init_pop_ab) // new seroconversion from infections
    - init_pop_ab * gamma // waning
    + delta * vacc[1];  // vaccination
  for (i in 2:t) {
    pop_ab[i] =
      pop_ab[i - 1]
      + (1 - beta) * infections[i] / ((1 - pop_ab[i - 1])) // new seroconversions from infections
      - pop_ab[i - 1] * gamma + // waning;
      delta * vacc[i]; // vaccination
  }
  return(pop_ab);
}
