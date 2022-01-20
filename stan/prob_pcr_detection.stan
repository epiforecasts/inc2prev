functions {
#include functions/detection_prob.stan
}
data {
#include chunks/data/pcr_testing.stan
}

parameters {
  real inc_mean;
  real<lower = 0> inc_sd;
  vector <lower = 0, upper = 1> [pcr_p] inf_at; 
  vector[3] pcr_eff;
  real<lower = 0> pcr_change;
}

model {
  // Priors on the incubation period
  inc_mean ~ normal(inc_mean_p[1], inc_mean_p[2]);
  inc_sd ~ normal(inc_sd_p[1], inc_sd_p[2]);

  // Prior on time infected at (as a proportion of time
  // from start to upper bound)
  inf_at ~ beta(3, 1);

  // Priors on piecewise linear (on logit) probability
  pcr_eff ~ std_normal();
  pcr_change ~ normal(5, 5) T[0, ];

  {
    // Infection time
    vector[pcr_n] tinf = pcr_inf_upper_bound .* inf_at;

    // Symptom onset likelihood
    // Probability of onset before first symptoms minus probability prior to
    // last asymptomatic test (or 0 if occurred prior to infection).
    vector[pcr_n] inf_to_sym = pcr_sym_at_test - tinf;
    vector[pcr_n] inf_to_lasym = pcr_last_asym_at_test - tinf;
    for (i in 1:pcr_n) {
      target += log(
        lognormal_cdf(inf_to_sym[i] | inc_mean, inc_sd) - 
        (inf_to_lasym[i] <= 0 ? 
          0 : lognormal_cdf(inf_to_lasym[i] | inc_mean, inc_sd)
        )
      );
    }

    // detection likelihood
    vector[pcr_p] inf_to_test = pcr_test_day - tinf[pcr_id];
    vector[pcr_p] pcr_p_d = detection_prob_logit(
      inf_to_test, pcr_eff, pcr_change
    );
    pcr_result ~ bernoulli_logit(pcr_p_d);
    // Add negative test at date of infection to constrain PCR detection
    target += bernoulli_logit_lupmf(0 | pcr_eff[1] - pcr_change * pcr_eff[2]);
  }
}

generated quantities {
  vector[301] pb;
  vector[301] k;
  for(j in 1:301) {
    k[j] = (j * 0.1) - 0.1;
  } 
  pb = detection_prob_logit(k, pcr_eff, pcr_change);
  pb = inv_logit(pb);
}
