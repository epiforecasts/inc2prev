library(data.table)

i2p_gp_tune_model <- function(path) {
  if (missing(path)) {
    path <- "stan/tune_inv_gamma.stan"
  }
  rstan::stan_model(path)
}

# define required stan data
i2p_data <- function(prev, ab, vacc, init_ab,
                     prob_detectable, unobserved_time = 14, horizon = 0,
                     inf_ab_delay = c(rep(0, 7 * 4), rep(1 / 7, 7)),
                     vacc_ab_delay = c(rep(0, 7 * 4), rep(1 / 7, 7)),
                     prop_seroconvert = c(2, 1), # 90%
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
                     ab_likelihood = TRUE,
                     var_col = NULL) {

  # check PMFS
  if (sum(inf_ab_delay) - 1 >= 1e-4) {
    stop("inf_ab_delay must sum to 1 rather than ", sum(inf_ab_delay))
  }

  if (sum(vacc_ab_delay) - 1 >= 1e-4) {
    stop("inf_ab_delay must sum to 1 rather than ", sum(vacc_ab_delay))
  }
  # extract a prevalence and build features
  prev <- data.table(prev)
  if (is.null(var_col)) {
    prev[, `..variable` := "dummy"]
  } else {
    setnames(prev, var_col, "..variable")
  }

  prev <- prev[, .(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    prev = middle,
    sd = (upper - lower) / (2 * 1.96),
    population,
    `..variable`
  )]
  prev[, `:=`(
    stime = as.integer(start_date - min(start_date)),
    etime = as.integer(end_date - min(start_date))
  )]
  data_prev <- t(as.matrix(
      dcast(prev, start_date ~ `..variable`, value.var = "prev")[, -1]
  ))
  data_prev_sd <- t(as.matrix(
      dcast(prev, start_date ~ `..variable`, value.var = "sd")[, -1]
  ))^2
  if (!is.null(ab)) {
    ## extract antibody prevalence and build features
    ab <- data.table(ab)
    if (is.null(var_col)) {
      ab[, `..variable` := "dummy"]
    } else {
      setnames(ab, var_col, "..variable")
    }
    ab <- ab[, .(
      start_date = as.Date(start_date),
      end_date = as.Date(end_date),
      prev = middle,
      sd = (upper - lower) / (2 * 1.96),
      `..variable`
    )]
    data_ab <- t(as.matrix(
      dcast(ab, start_date ~ `..variable`, value.var = "prev")[, -1]
    ))
    data_ab_sd <- t(as.matrix(
      dcast(ab, start_date ~ `..variable`, value.var = "sd")[, -1]
    ))^2
    ab_index <- match(rownames(data_ab), rownames(data_prev))
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
  vacc_base <- expand.grid(date = all_dates,
                          `..variable` = unique(prev$`..variable`))
  if (is.null(vacc)) {
    vacc <- data.table(vacc_base)[, vaccinated := 0]
  } else {
    vacc <- data.table(vacc)
    if (is.null(var_col)) {
      vacc[, `..variable` := "dummy"]
    } else {
      setnames(vacc, var_col, "..variable")
    }
    vacc <- vacc[, .(
      date = as.Date(date),
      vaccinated = vaccinated,
      ..variable
    )]
    setkey(vacc, date, `..variable`)
    vacc <- vacc[J(vacc_base), roll = 0]
    vacc <- vacc[is.na(vaccinated), vaccinated := 0]
  }
  data_vacc <- t(as.matrix(
    dcast(vacc, date ~ `..variable`, value.var = "vaccinated")[, -1]
  ))

  if (!is.null(init_ab)) {
    init_ab <- data.table(init_ab)
    if (is.null(var_col)) {
      init_ab[, `..variable` := "dummy"]
    } else {
      setnames(init_ab, var_col, "..variable")
    }
    init_ab <- init_ab[, .(
      prev = mean,
      sd = (upper - lower) / (2 * 1.96),
      `..variable`
    )]
  }

  ## summarise prob_detectable for simplicity
  prob_detectable <- melt(
    copy(prob_detectable),
    value.name = "p", id.vars = "sample"
  )
  prob_detectable[, time := as.numeric(as.character(variable))]
  prob_detectable <- prob_detectable[, .(
    median = median(p),
    mean = mean(p),
    sd = sd(p)
  ),
  by = time
  ]
  prob_detectable <- prob_detectable[,
    purrr::map(.SD, signif, digits = 3),
    .SDcols = c("mean", "median", "sd"),
    by = time
  ]
  # define baseline incidence
  init_inc_mean <- array(apply(data_prev, 1, mean) / sum(prob_detectable$mean))

  # build stan data
  dat <- list(
    ut = unobserved_time,
    t = length(all_dates),
    n = nrow(data_prev),
    n_ab = ifelse(is.null(ab), 0L, nrow(data_ab)),
    obs = ncol(data_prev),
    prev = data_prev,
    prev_sd2 = data_prev_sd,
    prev_stime = unique(prev$stime),
    prev_etime = unique(prev$etime),
    prob_detect_mean = rev(prob_detectable$mean),
    prob_detect_sd = rev(prob_detectable$sd),
    pbt = max(prob_detectable$time) + 1,
    init_inc_mean = logit(init_inc_mean),
    pbeta = prop_seroconvert,
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

  dat$ab_obs <- ifelse(is.null(ab), 0L, ncol(data_ab))
  if (!is.null(ab)) {
    dat <- c(dat, list(
      ab = data_ab,
      ab_sd2 = data_ab_sd^2,
      ab_stime = unique(ab$stime),
      ab_etime = unique(ab$etime),
      init_ab_mean = array(init_ab$prev),
      init_ab_sd = array(init_ab$sd),
      ab_index = array(ab_index),
      vacc = data_vacc
    ))
  } else {
    dat <- c(dat, list(
      ab = numeric(0),
      ab_sd2 = numeric(0),
      ab_stime = numeric(0),
      ab_etime = numeric(0),
      init_ab_mean = numeric(0),
      init_ab_sd = numeric(0),
      ab_index = integer(0),
      vacc = numeric(0)
    ))
  }

  # gaussian process parameters
  dat$M <- ceiling(dat$t * gp_m)
  dat$L <- 2
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$t
  }
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
      eta = array(
        rnorm(dat$M * dat$n, mean = 0, sd = 0.1), dim = c(dat$n, dat$M)
      ),
      alpha = array(
        truncnorm::rtruncnorm(dat$n, mean = 0, sd = 0.1, a = 0)
      ),
      rho = array(
        truncnorm::rtruncnorm(dat$n, mean = 36, sd = 21, a = 14, b = 90)
      ),
      sigma = array(truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)),
      prob_detect = purrr::map2_dbl(
        dat$prob_detect_mean, dat$prob_detect_sd / 10,
        ~ truncnorm::rtruncnorm(1, a = 0, b = 1, mean = .x, sd = .y)
      )
    )
    init_list$init_inc <- array(rnorm(dat$n, dat$init_inc_mean, 0.1))

    if (dat$diff_order > 0) {
      init_list$init_growth <- array(
        rnorm(dat$n * dat$diff_order, 0, 0.01), dim = c(dat$n, dat$diff_order)
      )
    }

    if (dat$ab_obs > 0) {
      init_list <- c(init_list, list(
        beta = array(inv_logit(rnorm(1, 2, 1))),
        gamma = array(inv_logit(rnorm(2, -9, 4))),
        delta = array(inv_logit(rnorm(1, 3, 1))),
        k = array(exp(rnorm(1, 0, 0.1))),
        l = array(exp(rnorm(1, 0, 0.1))),
        ab_sigma = array(
          truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)
        ),
        init_dab = array(truncnorm::rtruncnorm(
          dat$n, mean = dat$init_ab_mean, sd = dat$init_ab_sd / 10, a = 0
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
                      compile = TRUE, threads = FALSE, verbose = TRUE,
                      optimise = TRUE, ...) {
  if (missing(include)) {
    include <- "stan/functions"
  }
  if (optimise) {
    stanc_options <- list("O1")
  } else {
    stanc_options <- list()
  }

  if (compile) {
    if (verbose) {
      model <- cmdstanr::cmdstan_model(model,
        include_path = include,
        cpp_options = list(
          stan_threads = threads
        ),
        stanc_options = stanc_options,
        ...
      )
    } else {
      suppressMessages(
        model <- cmdstanr::cmdstan_model(model,
          include_path = include,
          cpp_options = list(
            stan_threads = threads
          ),
          stanc_options = stanc_options,
	  ...
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
