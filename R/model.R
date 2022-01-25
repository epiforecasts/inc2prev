library(data.table)

i2p_gp_tune_model <- function(path) {
  if (missing(path)) {
    path <- "stan/tune_inv_gamma.stan"
  }
  rstan::stan_model(path)
}

# define required stan data
i2p_data <- function(prev, ab, vacc, init_ab,
                     pb_params, unobserved_time = 40, horizon = 0,
                     init_cum_infections = c(0, 0),
                     inf_ab_delay = c(rep(0, 7 * 4), rep(1 / 7, 7)),
                     vacc_ab_delay = c(rep(0, 7 * 4), rep(1 / 7, 7)),
                     prop_dont_seroconvert = c(-2, 1), # 10%
                     inf_waning_rate = c(-9, 4),
                     vac_waning_rate = c(-9, 4), # 0.1%
                     vaccine_efficacy = c(3, 1), # 95%
                     gt = list(
                       mean = 3.64, mean_sd = 0.71, sd = 3.08,
                       sd_sd = 0.77, max = 15
                     ),
                     gp_m = 0.3, gp_ls = c(14, 90),
                     gp_tune_model = NULL,
                     differencing = 0,
                     prev_likelihood = TRUE,
                     ab_likelihood = TRUE) {

  # check PMFS
  if (sum(inf_ab_delay) - 1 >= 1e-4) {
    stop("inf_ab_delay must sum to 1 rather than ", sum(inf_ab_delay))
  }

  if (sum(vacc_ab_delay) - 1 >= 1e-4) {
    stop("inf_ab_delay must sum to 1 rather than ", sum(vacc_ab_delay))
  }
  # extract a prevalence and build features
  prev <- data.table(prev)[, .(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    prev = middle,
    sd = (upper - lower) / (2 * 1.96),
    population
  )]
  prev[, `:=`(
    stime = as.integer(start_date - min(start_date)),
    etime = as.integer(end_date - min(start_date))
  )]
  if (!is.null(ab)) {
    ## extract antibody prevalence and build features
    ab <- data.table(ab)[, .(
      start_date = as.Date(start_date),
      end_date = as.Date(end_date),
      prev = middle,
      sd = (upper - lower) / (2 * 1.96)
    )]
    model_start_date <- min(prev$start_date, ab$start_date)
    model_end_date <- max(prev$end_date, ab$end_date)
    ab[, `:=`(
      stime = as.integer(start_date - model_start_date),
      etime = as.integer(end_date - model_start_date)
    )]
  } else {
    model_start_date <- min(prev$start_date)
    model_end_date <- max(prev$end_date)
  }
  model_end_date <- model_end_date + days(horizon)
  all_dates <- seq(
    model_start_date - days(unobserved_time), model_end_date, by = "days"
  )
  prev[, `:=`(
    stime = as.integer(start_date - model_start_date),
    etime = as.integer(end_date - model_start_date)
  )]
  ## extract vaccination prevalence and fill missing dates with zeroes
  if (!is.null(vacc)) {
    vacc <- data.table(vacc)[, .(
      date = as.Date(date),
      vaccinated = vaccinated
    )]
    setkey(vacc, date)
    vacc <- vacc[J(all_dates), roll = 0]
    vacc <- vacc[is.na(vaccinated), vaccinated := 0]
  }
  if (!is.null(init_ab)) {
    init_ab <- data.table(init_ab)[, .(
      prev = mean,
      sd = (upper - lower) / (2 * 1.96)
    )]
  }

  # define baseline incidence
  init_inc_mean <- mean(prev$prev) / 
    (pb_params[variable %in% "mean_pb"]$mean * 
      pb_params[variable %in% "time"]$mean)

  # build stan data
  dat <- list(
    ut = unobserved_time,
    t = length(all_dates),
    obs = length(prev$prev),
    prev = prev$prev,
    prev_sd2 = prev$sd^2,
    prev_stime = prev$stime,
    prev_etime = prev$etime,
    pcr_eff_m = pb_params[grepl("beta", variable)]$mean,
    pcr_eff_sd = pb_params[grepl("beta", variable)]$sd,
    pcr_change_m = pb_params[variable %in% "cutpoint"]$mean,
    pcr_change_sd = pb_params[variable %in% "cutpoint"]$sd,
    prob_detect_mean = pb_params[grepl("prob_detect", variable)]$mean,
    prob_detect_sd = pb_params[grepl("prob_detect", variable)]$sd,
    pbt = pb_params[variable %in% "time"]$mean,
    pb_mode = 0,
    init_inc_mean = logit(init_inc_mean),
    pbeta = prop_dont_seroconvert,
    pgamma_mean = c(inf_waning_rate[1], vac_waning_rate[1]),
    pgamma_sd = c(inf_waning_rate[2], vac_waning_rate[2]),
    pdelta = vaccine_efficacy,
    linf_ab_delay = length(inf_ab_delay),
    inf_ab_delay = rev(inf_ab_delay),
    lvacc_ab_delay = length(vacc_ab_delay),
    vacc_ab_delay = rev(vacc_ab_delay),
    prev_likelihood = as.numeric(prev_likelihood),
    ab_likelihood = as.numeric(ab_likelihood)
  )

  dat$ab_obs <- ifelse(is.null(ab), 0L, length(ab$prev))
  if (!is.null(ab)) {
    dat <- c(dat, list(
      ab = ab$prev,
      ab_sd2 = ab$sd^2,
      ab_stime = ab$stime,
      ab_etime = ab$etime,
      init_ab_mean = init_ab$prev,
      init_ab_sd = init_ab$sd,
      vacc = vacc$vaccinated
    ))
  } else {
    dat <- c(dat, list(
      ab = numeric(0),
      ab_sd2 = numeric(0),
      ab_stime = numeric(0),
      ab_etime = numeric(0),
      init_ab_mean = numeric(0),
      init_ab_sd = numeric(0),
      vacc = numeric(0)
    ))
  }

  # gaussian process parameters
  dat$M <- ceiling(dat$t * gp_m)
  dat$L <- 2
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$t
  }
  gp_ls
  lsp <- tune_inv_gamma(gp_ls[1], gp_ls[2], gp_tune_model)
  dat$lengthscale_alpha <- lsp$alpha
  dat$lengthscale_beta <- lsp$beta
  dat$diff_order <- differencing

  # define generation time
  dat$gtm <- unlist(gt[c("mean", "mean_sd")])
  dat$gtsd <- unlist(gt[c("sd", "sd_sd")])
  dat$gtmax <- unlist(gt[c("max")])
  # nolint end
  return(dat)
}

library(truncnorm)
library(purrr)

i2p_inits <- function(dat) {
  inits <- function() {
    init_list <- list(
      eta = array(rnorm(dat$M, mean = 0, sd = 0.1)),
      alpha = array(truncnorm::rtruncnorm(1, mean = 0, sd = 0.1, a = 0)),
      sigma = array(truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)),
      rho = array(truncnorm::rtruncnorm(1, mean = 36, sd = 21, a = 14, b = 90)),
      pcr_eff = array(
        purrr::map2(dat$pcr_eff_m, dat$pcr_eff_sd, ~ rnorm(1, .x, .y * 0.1))
      ),
      pcr_change = rnorm(1, dat$pcr_change_m, dat$pcr_change_sd * 0.1),
      pb_sigma = abs(rnorm(1, 0.025, 0.001))
    )
    init_list$init_inc <- rnorm(1, dat$init_inc_mean, 0.1)

    if (dat$diff_order > 0) {
      init_list$init_growth <- array(rnorm(dat$diff_order, 0, 0.01))
    }

    if (dat$ab_obs > 0) {
      init_list <- c(init_list, list(
        beta = array(inv_logit(rnorm(1, -2, 0.4))),
        gamma = array(inv_logit(rnorm(2, -9, 0.4))),
        delta = array(inv_logit(rnorm(1, 3, 0.4))),
        ab_sigma = array(
          truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)
        ),
        init_dab = array(truncnorm::rtruncnorm(
          1, mean = dat$init_ab_mean, sd = dat$init_ab_sd / 10, a = 0
        ))
      ))
    }
    return(init_list)
  }
  return(inits)
}

#' Load and compile the model
#'
#' @param model A character string indicating the path to the model.
#' If not supplied the package default model is used.
#'
#' @param include A character string specifying the path to any stan
#' files to include in the model. If missing the package default is used.
#' @param compile Logical, defaults to `TRUE`. Should the model
#' be loaded and compiled using [cmdstanr::cmdstan_model()].
#'
#' @param threads Logical, defaults to `FALSE`. Should the model compile with
#' support for multi-thread support in chain. Note that this requires the use of
#' the `threads_per_chain` argument when model fitting using [enw_sample()],
#' and [epinowcast()].
#'
#' @param verbose Logical, defaults to `TRUE`. Should verbose
#' messages be shown.
#' @param ... Additional arguments passed to [cmdstanr::cmdstan_model()].
#'
#' @return A `cmdstanr` model.
#'
#' @family model
#' @export
#' @importFrom cmdstanr cmdstan_model
#' @examplesIf interactive()
#' mod <- i2p_model()
i2p_model <- function(model = "stan/inc2prev.stan", include,
                      compile = TRUE, threads = FALSE, verbose = TRUE, ...) {
  if (missing(include)) {
    include <- "stan/functions"
  }

  if (compile) {
    if (verbose) {
      model <- cmdstanr::cmdstan_model(model,
        include_path = include,
        cpp_options = list(
          stan_threads = threads
        ),
        ...
      )
    } else {
      suppressMessages(
        model <- cmdstanr::cmdstan_model(model,
          include_path = include,
          cpp_options = list(
            stan_threads = threads
          ), ...
        )
      )
    }
  }
  return(model)
}

#' Fit a CmdStan model using NUTS
#'
#' @param data A list of data as produced by [enw_as_data_list()].
#'
#' @param model A `cmdstanr` model object as loaded by [enw_model()].
#'
#' @param diagnostics Logical, defaults to `TRUE`. Should fitting diagnostics
#' be returned as a `data.frame`.
#'
#' @param ... Additional parameters passed to the `sample` method of `cmdstanr`.
#'
#' @return A `data.frame` containing the `cmdstanr` fit, the input data, the
#' fitting arguments, and optionally summary diagnostics.
#'
#' @family model
#' @export
#' @importFrom cmdstanr cmdstan_model
#' @importFrom posterior rhat
i2p_sample <- function(data, model = i2p_model(),
                       diagnostics = TRUE, ...) {
  fit <- model$sample(data = data, ...)

  out <- data.table(
    fit = list(fit),
    data = list(data),
    fit_args = list(list(...))
  )

  if (diagnostics) {
    diag <- fit$sampler_diagnostics(format = "df")
    diagnostics <- data.table(
      samples = nrow(diag),
      max_rhat = round(max(
        fit$summary(
          variables = NULL, posterior::rhat,
          .args = list(na.rm = TRUE)
        )$`posterior::rhat`,
        na.rm = TRUE
      ), 2),
      divergent_transitions = sum(diag$divergent__),
      per_divergent_transitions = sum(diag$divergent__) / nrow(diag),
      max_treedepth = max(diag$treedepth__)
    )
    diagnostics[, no_at_max_treedepth := sum(diag$treedepth__ == max_treedepth)]
    diagnostics[, per_at_max_treedepth := no_at_max_treedepth / nrow(diag)]
    out <- cbind(out, diagnostics)

    timing <- round(fit$time()$total, 1)
    out[, run_time := timing]
  }
  return(out[])
}
