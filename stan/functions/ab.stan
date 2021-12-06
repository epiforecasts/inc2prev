// Antibodies based on cumulative infections and vaccination
//
// @param infection Vector of new infections
//
// @param vacc Vector of new vaccinations
// 
// @param beta Proportion (0 - 1) that don't seroconvert
//
// @param gamma Daily rate of antibody waning (0 - 1). Leading to exponential
// waning.
//
// @param delta Vaccine efficacy (0 - 1). Assumes all or nothing antibodies.
//
// @param init_pop_ab Initial population that have antibodies (0 - 1)
//
// @param t Integer time index
//
// @examples
// source(here::here("R", "utils.R"))
// expose_stan_fns("ab.stan", "stan/functions")
//
// infections <- rep(10, 100)
// vacc <-  rep(1, 100)
// beta <- 0.05
// gamma <- 0.01
// delta <- 0.9
// init_pop_ab <- 0.2
// t <- 100
//
// detectable_antibodies(infections, vacc, beta, gamma, delta, init_pop_ab, t)
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
