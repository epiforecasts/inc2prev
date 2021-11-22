library("ggplot2")
library("here")
library("dplyr")
library("truncnorm")

## Get tools
source("R/utils.R")

geo_levels <- c("national", "regional", "local")
other_levels <- c("age_school")

levels <- c(geo_levels, other_levels)

prev <- read_cis()
estimates <- readRDS(here::here("outputs", "estimates.rds"))
samples <- readRDS(here::here("outputs", "samples.rds"))

early <- readRDS(here::here("data", "early.rds"))

for (level in levels) {
  ## 1) plot infections
  level_prev <- prev %>%
    filter(level == {{level}}) %>%
    mutate(variable = fct_inorder(variable))
  level_samples <- samples %>%
    filter(level == {{level}}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  level_estimates <- estimates %>%
    filter(level == {{level}}) %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  level_modelled_prev <- level_samples %>%
    filter(name == "pop_prev") %>%
    mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  if (level == "local") {
    level_early <- early %>%
      filter(level == "regional") %>%
      rename(region = variable) %>%
      inner_join(level_prev %>% select(-level), by = "region") %>%
      select(mean, low, high, variable)
  } else {
    level_early <- early %>% #
      filter(level == {{level}}) %>%
      mutate(variable = factor(variable, levels = levels(level_prev$variable)))
  }
  nvars <- n_distinct(level_prev$variable)
  p <- plot_prev(level_estimates, level_samples, level_prev)
  ggsave(here::here("figures", paste0("prev_", level, ".pdf")), p,
         width = 7 + 3 * floor(sqrt(nvars)),
         height = 7 + 3 * floor(sqrt(nvars)))
  p <- plot_trace(level_samples, "infections")
  ggsave(here::here("figures", paste0("inf_", level, ".pdf")), p,
         width = 7 + 3 * floor(sqrt(nvars)),
         height = 7 + 3 * floor(sqrt(nvars)))
  p <- plot_trace(level_samples, "cumulative_infections")
  ggsave(here::here("figures", paste0("cum_", level, ".pdf")), p,
         width = 7 + 3 * floor(sqrt(nvars)),
         height = 7 + 3 * floor(sqrt(nvars)))
  p <- plot_trace(level_samples, "R") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("R")
  ggsave(here::here("figures", paste0("R_", level, ".pdf")), p,
         width = 7 + 3 * floor(sqrt(nvars)),
         height = 7 + 3 * floor(sqrt(nvars)))
  p <- plot_trace(level_samples, "r")
  ggsave(here::here("figures", paste0("growth_", level, ".pdf")), p,
         width = 7 + 3 * floor(sqrt(nvars)),
         height = 7 + 3 * floor(sqrt(nvars)))
  early_samples <- level_early %>%
    group_by(variable) %>%
    summarise(rand = list(tibble(sample = 1:100,
                                 initial = rtruncnorm(n = 100, a = 0, mean = mean, sd = (high - low) / 4))),
              .groups = "drop") %>%
    unnest(rand)
  late_samples <- level_samples %>%
    filter(name == "cumulative_infections") %>%
    filter(date == max(date)) %>%
    pivot_longer(matches("[0-9]+"), names_to = "sample") %>%
    mutate(sample = as.integer(sample))
  combined <- late_samples %>%
    inner_join(early_samples, by = c("variable", "sample")) %>%
    mutate(value =  value + initial) %>%
    group_by(variable) %>%
    summarise(mean = mean(value),
              low = quantile(value, 0.025),
              high = quantile(value, 0.975),
              .groups = "drop")
  p <- ggplot(combined,
              aes(x = variable, y = mean, ymin = low, ymax = high)) +
    geom_point() +
    geom_linerange() +
    scale_y_continuous("Cumulative attack rate", labels = scales::percent_format(accuracy = 1L)) +
    theme_minimal() +
    expand_limits(y = 0) +
    xlab("")
  ggsave(here::here("figures", paste0("attack_rate_", level, ".pdf")),
         width = 10, height = 3)
}
