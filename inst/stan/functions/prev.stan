// Test using EpiNow2::expose_stan_fns("prev.stan", "stan/functions")
// convolve a pdf and case vector
vector detectable_cases(vector cases, vector[] pmf, int max_pmf, int t) {
    vector[t] conv_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
        conv_cases[s] += dot_product(cases[max(1, (s - max_pmf + 1)):s],
                                     tail(pmf[t], min(max_pmf, s)));
    }
   return(conv_cases);
  }
// Average observed cases across the period of time an estimate comes from
// observed_cases(1:100, c(1, 5, 6), c(4, 10, 9), 14, 3)
vector observed_cases(vector dcases, int[] prev_stime,
                      int[] prev_etime, int ut, int obs) {
    vector[obs] odcases = rep_vector(0, obs);

    for (i in 1:obs) {
        for (j in prev_stime[i]:prev_etime[i]) {
            odcases[i] += dcases[ut + j];
        }
        odcases[i] = odcases[i] / (prev_etime[i] - prev_stime[i]);
    }
    return(odcases);
}
