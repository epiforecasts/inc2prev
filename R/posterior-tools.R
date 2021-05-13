
#' Extract samples for a parameter with a date dimension from a Stan model
#'
#' @description Extracts a single parameter with a date dimension from a list
#' of `stan` output and returns it as a `data.table`.
#' @param fit An `rstanfit` object produced by `incidence()`
#' @param var Character string, variable to plot.
#' @param start_date A date, used to index the time plot.
#' @export
#' @return A data frame containing the parameter name, date, sample id and
#' sample value
extract_dated_parameter <- function(fit, var, start_date) {
  param_df <- as.data.table(
    t(
      as.data.table(
        extract(fit, pars = var)
      )
    )
  )
  param_df <- param_df[, time := 1:.N]
  param_df <- melt(param_df,
    id.vars = "time",
    variable.name = "var"
  )

  param_df <- param_df[, var := NULL][, sample := 1:.N, by = .(time)]
  param_df <- param_df[, date := start_date + time - 1, by = .(sample)]
  param_df <- param_df[, .(
    parameter = param, time, date,
    sample, value
  )]
  return(param_df)
}

#' Extract samples for a parameter with a date dimension from a Stan model
#'
#' @description Extracts a single parameter with a date dimension from a list
#' of `stan` output and returns it as a `data.table`.
#' @inheritParams extract_dated_parameter
#' @export
#' @return A data frame containing the parameter name, date, and summary
#' parameters
summarise_dated_parameter <- function(fit, var, start_date) {
  dt <- summary(fit, pars = var)$summary
  dt <- setDT(dt)
  dt <- dt[, time := 1:.N][date := start_date + time - 1]
  return(dt)
}


#' Extract samples from a parameter with a single dimension
#'
#' @inheritParams extract_dated_parameter
#' @export
#' @return A data frame containing the parameter name, sample id and sample
#' value.
extract_static_parameter <- function(param, samples) {
  data.table(
    parameter = param,
    sample = seq_length(samples[[param]]),
    value = samples[[param]]
  )
}
