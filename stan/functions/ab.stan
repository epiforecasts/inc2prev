// Antibodies based on cumulative infections and vaccination
//
// @param infection Vector of new infections
//
// @param vacc Vector of new vaccinations
// 
// @param beta Proportion (0 - 1) that seroconvert
//
// @param gamma Vector of daily rates of antibody waning (0 - 1) for both
// infection and vaccination. Assumes to exponential waning.
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
// infections <- rep(0.01, 100)
// vacc <-  rep(0.001, 100)
// beta <- 0.95
// gamma <- c(0.01, 0.001)
// delta <- 0.9
// init_pop_ab <- 0.2
// t <- 100
//
// detectable_antibodies(infections, vacc, beta, gamma, delta, init_pop_ab, t)
vector detectable_antibodies(vector infections, vector vacc,
                             real beta, vector gamma, real delta, real k, real l,
                             real init_pop_ab, int t) {
  vector[t] pop_ab;
  vector[t] inf_ab;
  vector[t] vac_ab;

  // Infection antibodies (assumes all previous from infection)
  inf_ab[1] = init_pop_ab // initial antibodies
    + beta * infections[1] // new seroconversion
    - init_pop_ab * gamma[1];
  // Vaccination antibodies
  vac_ab[1] = delta * vacc[1]; // vaccination
  // Population level antibodies
  pop_ab[1] = inf_ab[1] + vac_ab[1]; // population antibodies
  for (i in 2:t) {
    // Infection antibodies
    inf_ab[i] = inf_ab[i - 1] 
      + beta * infections[i] * pow(1 - pop_ab[i - 1], k) // new seroconversion
      - inf_ab[i - 1] * gamma[1]; // waning
    // Vaccination antibodies
    vac_ab[i] = vac_ab[i - 1] 
      + delta * vacc[i]  * pow(1 - pop_ab[i - 1], l) // new seroconversion
      - vac_ab[i - 1] * gamma[2]; // vaccination waning
    // Population antibodies
    pop_ab[i] = inf_ab[i] + vac_ab[i];
  }
  return(pop_ab);
}
