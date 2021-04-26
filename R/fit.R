library(truncnorm)
library(purrr)

stan_inits <- function(dat, n) {
  inits <- function() {
    list(
      eta = array(rnorm(dat$M, mean = 0, sd = 0.1)),
      alpha = array(truncnorm::rtruncnorm(1, mean = 0, sd = 0.1, a = 0)),
      sigma = array(truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)),
      rho = array(truncnorm::rtruncnorm(1, mean = 36, sd = 21, a = 14, b = 90))
    )
  }
  return(inits)
}

incidence <- function(dat, model, cores = 4, p, ...) {
  inits <- stan_inits(dat)

  fit <- model$sample(
    data = dat,
    init = inits,
    parallel_chains = cores,
    ...
  )
  if (!missing(p)) {
    p()
  }
  return(fit)
}

incidence_lapply <- function(dat_list, model, cores = 1, ...) {
  p <- progressr::progressor(along = dat_list)
  fits <- future_lapply(dat_list,
    incidence,
    model = model,
    cores = cores,
    p = p,
    future.seed = TRUE,
    ...
  )
  return(fits)
}

combine_incidence_fits <- function(incidence_list) {
  fits <- purrr::map(incidence_list, ~ read_stan_csv(.$output_files()))
  fit <- rstan::sflist2stanfit(fits)
  return(fit)
}
