#' Expose internal package stan functions in R
#'
#' @description This function exposes internal stan functions in R from a user
#' supplied list of target files. Allows for testing of stan functions in R and
#' potentially user use in R code.
#' @source https://github.com/epiforecasts/EpiNow2/
#' @param files A character vector indicating the target files
#' @param target_dir A character string indicating the target directory for the
#' file
#' @param ... Additional arguments passed to `rstan::expose_stan_functions`.
#' @return NULL
#' @export
#' @importFrom rstan expose_stan_functions stanc
#' @importFrom purrr map_chr
expose_stan_fns <- function(files, target_dir, ...) {
  functions <- paste0(
    "\n functions{ \n",
    paste(purrr::map_chr(
      files,
      ~ paste(readLines(file.path(target_dir, .)), collapse = "\n")
    ),
    collapse = "\n"
    ),
    "\n }"
  )
  expose_stan_functions(stanc(model_code = functions), ...)
  return(invisible(NULL))
}

#' Extract samples for a parameter with a date dimension from a Stan model
#'
#' @description Extracts a single parameter with a date dimension from a list
#' of `stan` output and returns it as a `data.table`.
#' @param fit An `rstanfit` object produced by `incidence()`
#' @param var Character string, variable to plot.
#' @param start_date A date, used to index the time plot.
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
#' @return A data frame containing the parameter name, sample id and sample value
extract_static_parameter <- function(param, samples) {
  data.table(
    parameter = param,
    sample = 1:length(samples[[param]]),
    value = samples[[param]]
  )
}
