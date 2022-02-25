// Code from: 
// https://github.com/epiforecasts/EpiNow2/tree/master/inst/stan/functions
// discretised truncated gamma pmf
vector discretised_gamma_pmf(array[] int y, real mu, real sigma, int max_val) {
  int n = num_elements(y);
  vector[n] pmf;
  real trunc_pmf;
  // calculate alpha and beta for gamma distribution
  real small = 1e-5;
  real large = 1e8;
  real c_sigma = sigma < small ? small : sigma;
  real c_mu = mu < small ? small : mu;
  real alpha = ((c_mu) / c_sigma)^2;
  real beta = (c_mu) / (c_sigma^2);
  // account for numerical issues
  alpha = alpha < small ? small : alpha;
  alpha = alpha > large ? large : alpha;
  beta = beta < small ? small : beta;
  beta = beta > large ? large : beta;
  // calculate pmf
  trunc_pmf = gamma_cdf(max_val + 1 | alpha, beta) - gamma_cdf(1 | alpha, beta);
  for (i in 1:n){
    pmf[i] = (gamma_cdf(y[i] + 1 | alpha, beta) - gamma_cdf(y[i] | alpha, beta)) /
    trunc_pmf;
  }
  return(pmf);
}
// calculate infectiousness (weighted sum of the generation time and infections)
// for a single time point
real update_infectiousness(vector infections, vector gt_pmf,
                           int seeding_time, int max_gt, int index){
  int inf_start = max(1, (index + seeding_time - max_gt));
  int inf_end = (index + seeding_time - 1);
  int pmf_accessed = min(max_gt, index + seeding_time - 1);
  real new_inf = dot_product(infections[inf_start:inf_end], tail(gt_pmf, pmf_accessed));
  return(new_inf);
}
// calculate Rt directly from inferred infections
vector calculate_Rt(vector infections, int seeding_time,
                    real gt_mean, real gt_sd, int max_gt,
                    int smooth) {
  vector[max_gt] gt_pmf;
  array[max_gt] int gt_indexes;
  int t = num_elements(infections);
  int ot = t - seeding_time;
  vector[ot] R;
  vector[ot] sR;
  vector[ot] infectiousness = rep_vector(1e-5, ot);
  // calculate PMF of the generation time
  for (i in 1:(max_gt)) {
    gt_indexes[i] = max_gt - i + 1;
  }
  gt_pmf = discretised_gamma_pmf(gt_indexes, gt_mean, gt_sd, max_gt);
  // calculate Rt using Cori et al. approach
  for (s in 1:ot) {
    infectiousness[s] += update_infectiousness(infections, gt_pmf, seeding_time,
                                               max_gt, s);
    R[s] = infections[s + seeding_time] / infectiousness[s];
  }
  if (smooth) {
    for (s in 1:ot) {
      real window = 0;
      sR[s] = 0;
      for (i in max(1, s - smooth):min(ot, s + smooth)) {
        sR[s] += R[i];
        window += 1;
      }
      sR[s] = sR[s] / window;
    }
  }else{
    sR = R;
  }
  return(sR);
}
