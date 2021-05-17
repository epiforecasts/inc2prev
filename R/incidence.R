library(data.table)

gt_opts <- function(mean, mean_sd, sd, sd_sd, max = 15) {
  dat <- list()
  dat$gtm <- c(mean, mean_sd)
  dat$gtsd <- c(sd, sd_sd)
  dat$gtmax <- max
  return(dat)
}

gp_opts <- function(m = 0.3, l = 2, ls_min = 9, ls_max = 90,
                    ls_tuner = inc2prev::load_model("tune_inv_gamma")){


}
#' Transform data into the format required for incidence modelling
#'
#' @description Transform available data into the format required for
#' `incidence()` and specify model options.
#'
#' @param observations
#' @param prob_detectable
#' @param ut Integer, defaults to 14 days. Initial unobserved period used to
#' seed the model prior to the first observed data point.
#' @param gt A list of options to specify the generation time. See `gt_opts()`
#' for detail
#' @param gp A list of options to specify the gaussian process. See `gp_opts()`
#' for details.
#' @param gp
#'
#' @return A function that generates a random sample of initial conditions
#' @export
#' @importFrom truncnorm rtruncnorm
#' @family incidence
incidence_data <- function(observations, prob_detectable, ut = 14,
                           population = 56286961,
                           observation_data = inc2prev::obs_data(),
                           additional_data = list(
                            gt = inc2prev::gt_opts(),
                            gp = inc2prev::gp_opts()
                          )) {
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

  # Combine additional options
  dat <- c(dat, additional_data)
  # gaussian process parameters
  dat$M <- ceiling(dat$t * gp$m)
  dat$L <- gp$l
  if (is.na(gp$ls_max) {
    gp$ls_max <- dat$t
  }

  if (missing(gp$ls_tuner)) {
    gp$ls_tuner <- load_model("tune_inv_gamma")
  }
  lsp <- tune_inv_gamma(gp$ls_tuner, gp$ls_min, gp$ls_max)
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


#' Initial conditions for the incidence model
#'
#' @description Produces random initial conditions used by stan
#' for the incidence estimation model.
#'
#' @param dat A list data as produced by `incidence_data()`. Must contain a
#' integer entry called M used to define the accuracy of the gaussian process
#' approximation.
#'
#' @return A function that generates a random sample of initial conditions
#' @export
#' @importFrom truncnorm rtruncnorm
#' @family incidence
#' @examples
#' dat <- list(M = 14)
#' inits <- incidence_inits(dat)
#' inits
#' inits()
incidence_inits <- function(dat) {
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


#' Incidence estimation model
#'
#' @description Fits an `inc2prev` model to recover incidence
#' from prevalence data using a probility of detection curve.
#' @param dat A list data as produced by `incidence_data()`.
#' @param model A stan model object as produced by `load_model("inc2prev")`
#' or a similar compiled stan model.
#' @param inits_fn A function that returns a function which produces
#' a set of initial conditions. Defaults to `incidence_inits()`.
#' @param fit_fn A function used to fit the `model` object. Defaults
#' to `rstan::sampling`.
#' @param p A `progressr` function used when fitting multiple models
#' to track progress
#' @param ... Additional arguments passed to `fit_fn`.
#' @return A fit stan model as produced by the option passed to
#' `fit_fn`.
#' @export
#' @family incidence
incidence <- function(dat, model = inc2prev::load_model("inc2prev"),
                      inits_fn = inc2prev::incidence_inits,
                      fit_fn = rstan::sampling, p, ...) {
  inits <- inits_fn(dat)

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

#' Fit multiple incidence models efficiently
#'
#' @description Fits multiple `incidence()` models in an efficient,
#' optionally in parallel, framework with user definable progression
#' settings.
#'
#' @details Progression information is displayed using the
#' `progressr` package. See the package documentation for options.
#' @param dat_list A list of lists of data as produced by `incidence_data()`.
#' @inheritParams incidence
#' @param ... Additional arguments passed to `incidence()` and/or
#' `future_lapply()`.
#' @return A fit stan model as produced by the option passed to
#' `fit_fn`.
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @family incidence
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

#' Combine the posteriors of multiple stan fit
#'
#' @description Combine the posteriors of multiple stan fits into
#' a single `rstan` `stanfit` object. This is likely useful when using
#' `incidence_lapply()` to fit models using multiple sampled probability
#' of detection curves. A simple wrapper around `rstan::sflist2stanfit()`
#' @details Progression information is displayed using the
#' `progresr` package. See the package documentation for options.
#' @param incidence_list A list of incidence model fits converted to be 
#' `stanfit` objects from th `rstan` packag
#' @return A `stanfit` model with combined posterior estimates
#' @export
#' @importFrom rstan sflist2stanfit
#' @family incidence
combine_incidence_fits <- function(incidence_list) {
  incidence <- sflist2stanfit(incidence_list)
  return(incidence)
}
