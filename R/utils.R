# plot trend with date over time

#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param fit DESCRIPTION.
#' @param var DESCRIPTION.
#' @param start_date DESCRIPTION.
#' @import data.table
#' @import ggplot2
#' @import data.table
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
plot_trend <- function(fit, var, start_date) {
  dt <- summary(fit, pars = var)$summary
  dt <- data.table::setDT(dt)
  dt <- dt[, time := 1:.N][date := start_date + time - 1]

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

plot_trace <- function(fit, var, start_date, samples = 100, alpha = 0.05,
                       rev_time = FALSE) {
  draws <- rstan::extract(fit)
  draws <- setDT(
    as.data.frame(draws[[var]])
  )

  as_draws_df() %>%
    as_tibble() %>%
    mutate(sample = 1:n()) %>%
    select(-.chain, -.iteration, -.draw) %>%
    pivot_longer(
      cols = -sample,
      names_to = "date",
      values_to = "value"
    ) %>%
    filter(sample <= samples) %>%
    group_by(sample) %>%
    mutate(time = 1:n(), date = date_start + time - 1)

  if (rev_time) {
    draws <- draws %>%
      mutate(date = rev(date))
  }
  draws <- draws %>%
    ungroup()

  plot <- draws %>%
    ggplot() +
    aes(x = date, y = value, group = sample) +
    geom_line(
      data = draws, alpha = alpha
    ) +
    theme_minimal() +
    labs(x = "Date")

  if (!rev_time) {
    plot <- plot +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %d")
  }
  return(plot)
}


library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

plot_prev <- function(fit, prev, samples = 100, date_start, alpha = 0.05,
                      data_source = "ONS Prevalence") {
  trace_plot <- plot_trace(
    fit,
    "pop_prev",
    date_start = date_start,
    samples = samples,
    alpha = alpha
  )

  summary_prev <- fit$summary(
    variables = "est_prev",
    ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))
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
