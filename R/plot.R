library(dplyr)
library(ggplot2)

# plot trend with date over time
plot_trend <- function(fit, var, date_start) {
  fit$summary(
    variables = var,
    ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))
  ) %>%
    mutate(
      time = 1:n(),
      date = date_start + time - 1
    ) %>%
    ggplot() +
    aes(x = date, y = `50%`, ymin = `5%`, ymax = `95%`) +
    geom_line(col = "lightblue", size = 1.4) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = "lightblue", size = 0.6
    ) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = NA, aes(ymin = `20%`, ymax = `80%`)
    ) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    theme_minimal()
}

plot_trace <- function(samples, var, alpha = 0.025) {
  long_samples <- samples %>%
    filter(name == var)

  plot <- long_samples %>%
    ggplot() +
    aes(x = date, y = value, group = sample) +
    geom_line(alpha = alpha) +
    theme_bw() +
    labs(x = "Date") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %d")
  return(plot)
}


library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

plot_prev <- function(estimates, samples, data, alpha = 0.05,
                      modelled = "dcases", observed = "est_prev",
                      data_source = "ONS Prevalence") {
  trace_plot <- plot_trace(
    samples,
    modelled,
    alpha = alpha
  )

  summary_prev <- estimates %>%
    filter(name == {{ observed }}) %>%
    mutate(
      middle = `50%`,
      lower = `5%`,
      upper = `95%`,
      type = "Modelled"
    ) %>%
    bind_rows(data %>%
      mutate(
        type = "Estimate"
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
      size = 1, alpha = 0.2
    ) +
    geom_point(
      data = summary_prev,
      aes(
        y = middle, ymin = NULL, ymax = NULL, group = NULL,
        col = type
      ), size = 1.1, alpha = 0.2
    ) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Dark2") +
    guides(col = guide_legend(title = data_source))
}

plot_ltla <- function(estimates, areas, names = c(), var = "pop_prev",
                      days = 60, var_name = "Prevalence") {
  estimates <- estimates %>%
    filter(name == {{ var }}) %>%
    filter(date > max(date) - {{ days }})
  if (length(names) > 0) {
    search_str <- paste0(
      "(",
      paste(names, collapse = "|"),
      ")"
    )
    areas <- areas %>%
      mutate(highlighted = grepl(search_str, ltla_name)) %>%
      group_by(geography_code) %>%
      summarise(highlighted = any(highlighted), .groups = "drop") %>%
      mutate(
        highlighted = if_else(highlighted, "yes", "no"),
        highlighted = factor(highlighted, levels = c("yes", "no"))
      )
  }
  estimates <- estimates %>%
    inner_join(areas %>% rename(variable = geography_code), by = "variable")
  aesthetics <- list(
    x = "date",
    y = "`50%`",
    group = "variable"
  )
  if (length(names) > 0) {
    aesthetics[["colour"]] <- "highlighted"
    aesthetics[["alpha"]] <- "highlighted"
  } else {
    aesthetics[["colour"]] <- "region"
  }
  p <- ggplot(estimates, mapping = do.call(aes_string, aesthetics)) +
    geom_line()
  if (length(names) > 0) {
    p <- p +
      scale_colour_manual("", values = c("red", "black")) +
      scale_alpha_manual("", values = c(1, 0.25)) +
      theme(legend.position = "none")
  } else {
    p <- p +
      scale_colour_brewer("Region", palette = "Paired")
  }
  p <- p +
    theme_minimal() +
    xlab("") +
    ylab(var_name)
  return(p)
}
