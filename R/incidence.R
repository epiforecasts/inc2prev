library(data.table)

# define required stan data
incidence_data <- function(observations, prob_detectable, ut = 14,
                           population = 56286961,
                           gt = list(
                             mean = 3.64, mean_sd = 0.71, sd = 3.08,
                             sd_sd = 0.77, max = 15
                           ),
                           gp_m = 0.3,
                           gp_l = 2,
                           gp_ls = c(14, 90),
                           gp_ls_tuner) {
  # nolint start
  # extract a single region for prevalence and build features
  prev <- copy(observations)
  prev <- prev[, .(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    date = as.Date(date),
    prev = middle,
    sd = (upper - lower) / (2 * 1.96)
  )]
  prev[, `:=`(
    time = as.integer(date - min(start_date)),
    stime = as.integer(start_date - min(start_date)),
    etime = as.integer(end_date - min(start_date))
  )]

  # boostrapped probability of detection
  pd_samples <- copy(prob_detectable)
  pd_samples <- split(pd_samples, by = "sample")
  pd_samples <- map(pd_samples, ~ .[, sample := NULL])

  # define common data
  dat <- list(
    ut = ut,
    ot = max(prev$etime),
    t = ut + max(prev$etime),
    obs = length(prev$prev),
    prev = prev$prev,
    prev_sd2 = prev$sd^2,
    prev_time = prev$time,
    prev_stime = prev$stime,
    prev_etime = prev$etime,
    N = population
  )

  # gaussian process parameters
  dat$M <- ceiling(dat$t * gp_m)
  dat$L <- gp_l
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$t
  }

  if (missing(gp_ls_tuner)) {
    gp_ls_tuner <- load_model("tune_inv_gamma")
  }
  lsp <- tune_inv_gamma(gp_ls_tuner, gp_ls[1], gp_ls[2])
  dat$lengthscale_alpha <- lsp$alpha
  dat$lengthscale_beta <- lsp$beta

  # define generation time
  dat$gtm <- unlist(gt[c("mean", "mean_sd")])
  dat$gtsd <- unlist(gt[c("sd", "sd_sd")])
  dat$gtmax <- unlist(gt[c("max")])

  dat_list <- map(pd_samples, function(pd) {
    if (is.null(pd$time)) {
      prob_detect <- rep(list(rev(unlist(pd))), dat$t)
    } else {
      prob_detect <- split(pd, by = "time")
      prob_detect <- map(pd_samples, ~ rev(unlist(.[, sample := NULL])))
    }
    # define baseline incidence
    baseline_inc <- prev$prev[1] * prob_detect[[1]][ut]

    # build stan data
    sample_dat <- list(
      prob_detect = prob_detect,
      pbt = length(prob_detect[[1]]),
      inc_zero = log(baseline_inc / (baseline_inc + 1))
    )

    return(c(dat, sample_dat))
  })

  if (length(dat_list) == 1) {
    return(dat_list[[1]])
  } else {
    return(dat_list)
  }
  # nolint end
}

library(truncnorm)
library(purrr)

incidence_inits <- function(dat, n) {
  inits <- function() {
    list(
      eta = array(rnorm(dat$M, mean = 0, sd = 0.1)),
      alpha = rtruncnorm(1, mean = 0, sd = 0.1, a = 0),
      sigma = rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0),
      rho = rtruncnorm(1, mean = 36, sd = 21, a = 14, b = 90)
    )
  }
  return(inits)
}



incidence <- function(dat, model,
                      inits_fn = incidence_inits,
                      fit_fn = rstan::sampling, p, ...) {
  inits <- inits_fn(dat)

  if (missing(model)) {
    model <- load_model("tune_inv_gamma")
  }

  fit <- do.call(fit_fn, list(
    object = model,
    data = dat,
    init = inits,
    ...
  ))
  if (!missing(p)) {
    p()
  }
  return(fit)
}

incidence_lapply <- function(dat_list, model, cores = 1, ...) {
  p <- progressor(along = dat_list)
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
  incidence <- sflist2stanfit(incidence_list)
  return(incidence)
}
