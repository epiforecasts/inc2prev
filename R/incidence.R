#' Incidence estimation model
#'
#' @description Fits an `inc2prev` model to recover incidence
#' from prevalence data using a probility of detection curve.
#' @param prev Observed positiviy prevalence data
#' @param ab Observed antibody prevalence data
#' @param vacc Observed vaccination data
#' @param init_ab Observed initial antibody data
#' @param prob_detection A dataframe of posterior samples for the probability
#' of detection
#' @param data_args A list of arguments to pass to `i2p_data()`
#' @param model A stan model object as produced by `i2p_model()`
#' or a similar compiled stan model.
#' @param inits_fn A function that returns a function which produces
#' a set of initial conditions. Defaults to `i2p_inits()`.
#' @param fit_fn A function used to fit the `model` object. Defaults
#' to `i2p_sample`.
#' @param p A `progressr` function used when fitting multiple models
#' to track progress
#' @param ... Additional arguments passed to `fit_fn`.
#' @return A fit stan model as produced by the option passed to
#' `fit_fn`.
#' @export
#' @family incidence
incidence <- function(prev, ab = NULL, vacc = NULL, init_ab = NULL, prob_detect,
                      data_args = list(),
                      model = i2p_model(),
                      variables = c(
                        "est_prev", "pop_prev", "infections",
                        "dcases", "r", "R"
                      ),
                      quantiles = seq(0.05, 0.95, by = 0.05),
                      samples = 100,
                      keep_fit = FALSE,
                      p, ...) {
  dat <- do.call(
    i2p_data, c(
      list(
        prev = prev,
        ab = ab,
        vacc = vacc,
        init_ab = init_ab,
        prob_detectable = prob_detect
      ),
      data_args
    )
  )
  inits <- i2p_inits(dat)

  fit <- do.call(i2p_sample, list(
    model = model,
    data = dat,
    init = inits,
    ...
  ))

  fit[, summary := list(
    i2p_summarise(fit[[1]], variables = variables, quantiles = quantiles) |>
      i2p_add_date(prev = prev, data = dat)
  )]

  fit[, samples := list(
    i2p_draws(fit[[1]], variables = variables, samples = samples) |>
      i2p_add_date(prev = prev, data = dat)
  )]
  if (!keep_fit) {
    fit[, fit := NULL]
  }

  if (!missing(p)) {
    p()
  }
  return(fit[])
}

#' Fit multiple incidence models efficiently
#'
#' @description Fits multiple `incidence()` models in an efficient,
#' optionally in parallel, framework with user definable progression
#' settings.
#'
#' @details Progression information is displayed using the
#' `progressr` package. See the package documentation for options.
#' @param prev_list A list of lists of data as produced by `i2p_data()`.
#' @inheritParams incidence
#' @param ... Additional arguments passed to `incidence()` and/or
#' `future_lapply()`.
#' @return A fit stan model as produced by the option passed to
#' `fit_fn`.
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @family incidence
incidence_lapply <- function(prev_list, ...) {
  p <- progressor(along = prev_list)
  fits <- future_lapply(prev_list,
    incidence,
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
