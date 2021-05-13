#' Plot method for inc2prev
#'
#' @description Plot method for models fit using the `inc2prev` package.
#' @param x A model fit using the `incidence()` method.
#' @param type A character vector indicating the plot method to use. Supported
#' options are "cri" (plot credible intervals), "trace" (plot posterior
#' samples), and "obs" (plot observations against the posterior).
#' @param ... Pass additional arguments to underlying plot functions
#' @aliases plot
#' @method plot inc2prev
#' @return A `ggplot2` object
#' @export
plot.inc2prev <- function(type = "trace", ...) {
  type <- match.arg(type, choices = c("cri", "trace", "obs"))

  do.call(paste0("plot_", type), ...)
}

#' Plot posterior credible intervals over time
#'
#' @inheritParams summarise_dated_parameter
#' @return A `ggplot2` object
plot_cri <- function(fit, var, start_date) {
  dt <- summarise_dated_parameter(
    fit, var, start_date
  )

  ggplot(dt) +
    aes(x = date, y = `50%`, ymin = `2.5%`, ymax = `97.5%`) +
    geom_line(col = "lightblue", size = 1.4) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = "lightblue", size = 0.6
    ) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = NA, aes(ymin = `25%`, ymax = `75%`)
    ) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    theme_minimal()
}


#' Plot posterior samples over time
#'
#' @param samples Integer, defaults to 100. The number of posterior
#' samples to extract.
#' @param alpha Numeric, defaults to 0.05. The alpha value to use in
#' the`ggplot2`  plot for the individual posterior traces.
#' @inheritParams extract_dated_parameter
#' @return A `ggplot2` object
plot_trace <- function(fit, var, start_date, samples = 100, alpha = 0.05) {
  draws <- extract_dated_parameter(
    fit = fit,
    var = var,
    start_date = start_date,
    samples = samples
  )

  plot <- draws %>%
    ggplot() +
    aes(x = date, y = value, group = sample) +
    geom_line(
      data = draws, alpha = alpha
    ) +
    theme_minimal() +
    labs(x = "Date")
  return(plot)
}

#' Plot observations against posterior samples over time
#'
#' @param data_source A character string indicating the label to use to
#' identify the observations.
#' @inheritParams incidence_data
#' @inheritParams plot_trace
#' @return A `ggplot2` object
#' @importFrom scales percent
plot_obs <- function(fit, observations, samples = 100, start_date,
                     alpha = 0.05, data_source) {
  trace_plot <- plot_trace(
    fit,
    "pop_prev",
    start_date = start_date,
    samples = samples,
    alpha = alpha
  )

  summary_prev <- summarise_dated_parameter(
    var = "est_prev",
    start_date = start_date,
  )

  summary_prev <- summary_prev[
    ,
    `:=`(
      middle = `50%`,
      lower = `5%`,
      upper = `95%`,
      date = observations$date,
      type = "Modelled"
    )
  ]
  summary_pre <- rbind(
    summary_prev,
    copy(observations)[, `:=`(type = "Estimate", date = prev$date)]
  )

  trace_plot +
    scale_y_continuous(labels = percent) +
    labs(y = "Prevalence", x = "Date") +
    geom_linerange(
      data = summary_prev,
      aes(
        y = NULL, ymin = lower, ymax = upper, group = NULL,
        col = type
      ),
      size = 1, alpha = 0.8
    ) +
    geom_point(
      data = summary_prev,
      aes(
        y = middle, ymin = NULL, ymax = NULL, group = NULL,
        col = type
      ), size = 1.1
    ) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Dark2") +
    guides(col = guide_legend(title = data_source))
}
