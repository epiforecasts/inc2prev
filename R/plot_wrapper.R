plot_wrapper <- function(level, prev, ab = NULL, samples, estimates, early,
                         dont_seroconvert = 0) {
  level_prev <- prev %>%
    filter(level == {{ level }}) %>%
    mutate(variable = fct_inorder(variable))
  if (nrow(level_prev) == 0) {
    stop("No data available with these filter settings")
  }
  if (!is.null(ab)) {
    level_ab <- ab %>%
      filter(level == {{ level }}) %>%
      mutate(variable = fct_inorder(variable))
    if (nrow(level_ab) == 0) {
      stop("No antibody data available with these filter settings")
    }
    suffix <- "_ab"
  } else {
    suffix <- ""
  }
  level_samples <- samples %>%
    filter(level == {{ level }}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  if (nrow(level_samples) == 0) {
    stop("No samples available with these filter settings")
  }
  level_estimates <- estimates %>%
    filter(level == {{ level }}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  if (nrow(level_estimates) == 0) {
    stop("No estimates available with these filter settings")
  }
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
  p <- plot_prev(level_estimates, level_samples, level_prev) +
    facet_wrap(~variable, scales = "free_y")
  ggsave(here::here("figures", paste0("prev_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  if (!is.null(ab)) {
    p <- plot_prev(level_estimates, level_samples, level_ab,
      modelled = "dab", observed = "est_ab"
    ) +
      facet_wrap(~variable)
    ggsave(here::here("figures", paste0("ab_", level, suffix, ".png")), p,
      width = 7 + 3 * floor(sqrt(nvars)),
      height = 2 + 3 * floor(sqrt(nvars))
    )
  }
  ## 2) plot incidence
  p <- plot_trace(level_samples, "infections") +
    facet_wrap(~variable, scales = "free_y")
  ggsave(here::here("figures", paste0("inf_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 3) plot R
  p <- plot_trace(level_samples, "R") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("R") +
    facet_wrap(~variable)
  ggsave(here::here("figures", paste0("R_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 4) plot growth
  p <- plot_trace(level_samples, "r") +
    facet_wrap(~variable)
  ggsave(here::here("figures", paste0("growth_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 5) plot cumulative incidence
  ## sample from late April seroprevalence
  if (!("cumulative_infections" %in% unique(level_samples$name))) {
    level_samples <- level_samples %>%
      filter(name == "infections") %>%
      arrange(sample, index, variable) %>%
      group_by(sample, variable) %>%
      mutate(value = cumsum(value)) %>%
      ungroup() %>%
      mutate(name = "cumulative_infections") %>%
      bind_rows(level_samples)
  }
  p <- plot_trace(level_samples, "cumulative_infections") +
    scale_y_continuous("Cumulative incidence",
      labels = scales::percent_format(1L)
    ) +
    xlab("") +
    facet_wrap(~variable)
  ggsave(here::here("figures", paste0("car_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 6) cumulative infections
  ## get samples of estimated final cumulative incidence
  cumulative_infection_samples <- level_samples %>%
    filter(name == "cumulative_infections")
  ## combine early seroprevalence with estimates of attack rates since start of CIS # nolint
  n_samples <- max(cumulative_infection_samples$sample)
  early_samples <- level_early %>%
    group_by(variable) %>%
    summarise(
      rand = list(tibble(
        sample = seq_len(n_samples),
        initial = rtruncnorm(
          n = n_samples, a = 0, mean = mean,
          sd = (upper - lower) / 4
        )
      )),
      .groups = "drop"
    ) %>%
    unnest(rand) %>%
    mutate(initial = initial / (1 - dont_seroconvert))
  combined_samples <- cumulative_infection_samples %>%
    inner_join(early_samples, by = c("variable", "sample")) %>%
    mutate(value = value + initial) %>%
    select(-initial)
  p <- plot_trace(combined_samples, "cumulative_infections") +
    scale_y_continuous("Cumulative incidence",
      labels = scales::percent_format(1L)
    ) +
    xlab("") +
    facet_wrap(~variable)
  ggsave(here::here("figures", paste0("cum_inc_", level, suffix, ".png")), p,
    width = 7 + 3 * floor(sqrt(nvars)),
    height = 2 + 3 * floor(sqrt(nvars))
  )
  ## 7) final attack rates
  ## combine early seroprevalence with estimates of attack rates
  ## since start of CIS
  combined_aggregate <- combined_samples %>%
    filter(date == max(date)) %>%
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
    theme_light() +
    expand_limits(y = 0) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave(here::here("figures", paste0("car_", level, suffix, ".png")),
    width = 7, height = 4
  )
  ## 8) final antibodies
  if (!is.null(ab)) {
    antibodies <- level_samples %>%
      filter(
        date == max(date),
        name == "dab"
      ) %>%
      group_by(variable) %>%
      summarise(
        mean = mean(value),
        low = quantile(value, 0.025),
        high = quantile(value, 0.975),
        .groups = "drop"
      )
    p <- ggplot(
      antibodies,
      aes(x = variable, y = mean, ymin = low, ymax = high)
    ) +
      geom_point() +
      geom_linerange() +
      scale_y_continuous(
        "Seroprevalence",
        labels = scales::percent_format(accuracy = 1L)
      ) +
      theme_light() +
      expand_limits(y = 0) +
      xlab("") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(here::here("figures", paste0("seroprevalence_", level, suffix, ".png")),
      width = 7, height = 4
    )
  }
  return(invisible(NULL))
}
