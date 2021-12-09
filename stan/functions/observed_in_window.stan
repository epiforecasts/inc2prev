// Average observed cases across the period of time an estimate comes from
// observed_in_window(1:100, c(1, 5, 6), c(4, 10, 9), 14, 3)
vector observed_in_window(vector dcases, int[] prev_stime,
                          int[] prev_etime, int ut, int obs) {
    vector[obs] odcases = rep_vector(0, obs);

    for (i in 1:obs) {
        for (j in prev_stime[i]:prev_etime[i]) {
            odcases[i] += dcases[ut + j];
        }
        odcases[i] = odcases[i] / (prev_etime[i] - prev_stime[i] + 1);
    }
    return(odcases);
}
