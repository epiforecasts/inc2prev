// Code from: 
// https://github.com/epiforecasts/EpiNow2/tree/master/inst/stan/functions
// Calculate growth rate
vector calculate_growth(vector infections, int seeding_time) {
  int t = num_elements(infections);
  int ot = t - seeding_time;
  vector[t] log_inf = log(infections); 
  vector[ot] growth = log_inf[(seeding_time + 1):t] - log_inf[seeding_time:(t - 1)];
  return(growth);
}
