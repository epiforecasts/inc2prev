plot_wrapper <- function(level, prev, samples, estimates, early) {
  level_prev <- prev %>%
    filter(level == {{ level }}) %>%
    mutate(variable = fct_inorder(variable))
  level_samples <- samples %>%
    filter(level == {{ level }}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  level_estimates <- estimates %>%
    filter(level == {{ level }}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  level_modelled_prev <- level_samples %>%
    filter(name == "pop_prev") %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  if (level == "local") {
    level_early <- early %>%
      filter(level == "regional") %>%
      rename(region = variable) %>%
      inner_join(level_prev %>%
        select(-level),
      by = "region"
      ) %>%
      select(mean, low, high, variable)
  } else {
    level_early <- early %>% #
      filter(level == {{ level }}) %>%
      mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  }
  nvars <- n_distinct(level_prev$variable)
  ## 1) plot prevalence
  p <- plot_prev(level_estimates, level_samples, level_prev)
  ggsave(here::here("figures", paste0("prev_", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 2) plot incidence
  p <- plot_trace(level_samples, "infections")
  ggsave(here::here("figures", paste0("inf_", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 3) plot R
  p <- plot_trace(level_samples, "R") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("R")
  ggsave(here::here("figures", paste0("R_", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 4) plot growth
  p <- plot_trace(level_samples, "r")
  ggsave(here::here("figures", paste0("growth_", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 5) plot cumulative incidence
  ## sample from late April seroprevalence
  early_samples <- level_early %>%
    group_by(variable) %>%
    summarise(
      rand = list(tibble(
        sample = 1:100,
        initial = rtruncnorm(n = 100, a = 0, mean = mean, sd = (high - low) / 4)
      )),
      .groups = "drop"
    ) %>%
    unnest(rand)
  ## get samples of estimated final cumulative incidence
  cumulative_infection_samples <- level_samples %>%
    filter(name == "cumulative_infections") %>%
    pivot_longer(matches("[0-9]+"), names_to = "sample") %>%
    mutate(sample = as.integer(sample))
  combined_samples <- cumulative_infection_samples %>%
    inner_join(early_samples, by = c("variable", "sample")) %>%
    mutate(value = value + initial) %>%
    select(-initial) %>%
    pivot_wider(names_from = "sample")
  p <- plot_trace(combined_samples, "cumulative_infections") +
    scale_y_continuous("Cumulative incidence",
      labels = scales::percent_format(1L)
    ) +
    xlab("")
  ggsave(here::here("figures", paste0("car_", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 6) plot final attack rate
  ## combine early seroprevalence with estimates of attack rates since start of CIS # nolint
  combined_aggregate <- combined_samples %>%
    filter(date == max(date)) %>%
    pivot_longer(matches("[0-9]+"), names_to = "sample") %>%
    group_by(variable) %>%
    summarise(
      mean = mean(value),
      low = quantile(value, 0.025),
      high = quantile(value, 0.975),
      .groups = "drop"
    )
  p <- ggplot(
    combined_aggregate,
    aes(x = variable, y = mean, ymin = low, ymax = high)
  ) +
    geom_point() +
    geom_linerange() +
    scale_y_continuous(
      "Cumulative attack rate",
      labels = scales::percent_format(accuracy = 1L)
    ) +
    theme_minimal() +
    expand_limits(y = 0) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave(here::here("figures", paste0("attack_rate_", level, ".png")),
    width = 7, height = 4
  )
  ## scatter attack rate against current incidence
  p <- level_estimates %>%
    filter(name == "infections") %>%
    filter(date == max(date) - 7) %>%
    select(variable, median_inc = `50%`) %>%
    left_join(combined_aggregate %>%
      select(variable, mean_car = mean), by = "variable") %>%
    left_join(level_prev %>%
      select(variable, population) %>%
      distinct(), by = "variable") %>%
    mutate(median_inc = median_inc / population) %>%
    ggplot(aes(x = mean_car, median_inc)) +
    geom_jitter() +
    ylab("Latest incidence") +
    xlab("Cumulative attack rate") +
    theme_minimal()
  ggsave(here::here("figures", paste0("car_vs_incidence", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 7 + 3 * floor(sqrt(nvars))
  )
  p <- level_estimates %>%
    filter(name == "R") %>%
    filter(date == max(date) - 7) %>%
    select(variable, median_R = `50%`) %>%
    left_join(combined_aggregate %>%
      select(variable, mean_car = mean), by = "variable") %>%
    ggplot(aes(x = mean_car, median_R)) +
    geom_jitter() +
    ylab("Latest reproduction number") +
    xlab("Cumulative attack rate") +
    theme_minimal()
  ggsave(here::here("figures", paste0("car_vs_R", level, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 7 + 3 * floor(sqrt(nvars))
  )
  return(invisible(NULL))
}
