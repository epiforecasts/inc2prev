
#' Plot summarised trend over time in posteriors
#'
#' @inheritparams extract_parameter
#' @return A `ggplot2` object
plot_trend <- function(fit, var, start_date) {
  dt <- summarise_dated_parameter(fit, var, start_date)

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
#' @inheritParams extracted_dated_parameter
#' @return A `ggplot2` object
#' @importFrom rstan extract
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


plot_prev <- function(fit, prev, samples = 100, date_start, alpha = 0.05,
                      data_source = "ONS Prevalence") {
  trace_plot <- plot_trace(
    fit,
    "pop_prev",
    date_start = date_start,
    samples = samples,
    alpha = alpha
  )

  summary_prev <- summarise_dated_parameter(
    var = "est_prev",
    start_date = prev_date + 2
  ) %>%
    mutate(
      middle = `50%`,
      lower = `5%`,
      upper = `95%`,
      date = prev$date + 2,
      type = "Modelled"
    ) %>%
    bind_rows(prev %>%
      mutate(
        type = "Estimate",
        date = date - 2
      ))

  trace_plot +
    scale_y_continuous(labels = scales::percent) +
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
