// convolve a pdf and case vector
vector detectable_cases(vector cases, vector pmf, int max_pmf, int t) {
    vector[t] conv_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
        conv_cases[s] += dot_product(cases[max(1, (s - max_pmf + 1)):s],
                                     tail(pmf, min(max_pmf, s)));
    }
   return(conv_cases);
  }

vector observed_cases(vector dcases, int[] prev_time, int ut, int obs) {
    vector[obs] odcases;

    for (i in 1:obs) {
        odcases[i] = dcases[ut + prev_time[i]];
    }
    return(odcases);
}
