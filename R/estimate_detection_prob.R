detection_prob_obs_to_list <- function(dt) {

if (!is.null(dt)) {
  ids <- unique(dt[,
   .(num_id, first_sym_day, last_asym_day, inf_upper_bound)]
  )

  dat <- list(
    pb_p = nrow(ids),
    pb_n = nrow(dt),
    pb_id = dt$num_id,
    pb_test_day = dt$day,
    pb_result = as.numeric(dt$result),
    pb_sym_at_test = ids$first_sym_day,
    pb_last_asym_at_test = ids$last_asym_day,
    pb_inf_upper_bound = ids$inf_upper_bound
  )
}else{
  dat <- list(
    pb_n = 0,
    pb_p = 0,
    pb_id = numeric(0),
    pb_test_day = numeric(0),
    pb_result = numeric(0),
    pb_sym_at_test = numeric(0),
    pb_last_asym_at_test = numeric(0),
    pb_inf_upper_bound = numeric(0)
  )
}
  return(dat)
}

get_inc_period <- function(inc_mean = c(1.621, 0.0640),
                           inc_sd = c(0.418, 0.0691)) {
  list(
    inc_mean_p = inc_mean,
    inc_sd_p = inc_sd
  )
}

estimate_detection_prob <- function(obs, inc_period = get_inc_period(), ...) {
  dt <- detection_prob_obs_to_list(obs)
  dt <- c(dt, inc_period)

  mod <- cmdstanr::cmdstan_model(
    here::here("stan", "prob_detection.stan"),
    include_paths = here::here("stan")
  )
  fit <- mod$sample(
    data = dt,
    ...
  )
  return(fit)
}