library(here)
library(purrr)
library(ggplot2)
library(socialmixr)
library(dplyr)
library(readr)
library(scales)
library(tidyr)

source(here::here("R", "read.R"))
prev <- read_cis()
local_region <- prev %>%
  filter(level == "local") %>%
  select(level, variable, region) %>%
  distinct()

dir.create(here::here("figures", "additional"),
           showWarnings = FALSE, recursive = TRUE)

## Get tools
files <-
  list.files(here("outputs"), pattern = "estimates_[^_]+\\.csv", full.names = TRUE)
data <- map(files, read_csv) %>%
  bind_rows() %>%
  separate(variable, c("variable", "geography"), sep = "\\|")
data <- left_join(data, local_region, by = c("level", "variable"))

cum_files <-
  list.files(here("outputs"), pattern = "cumulative_[^_]+\\.rds", full.names = TRUE)
if (length(cum_files) > 0) {
  cum_data <- map(cum_files, readRDS) %>%
    bind_rows() %>%
    group_by(name, date, variable, variant, level) %>%
    summarise(
      q5 = quantile(value, 0.05),
      q95 = quantile(value, 0.95),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      q45 = quantile(value, 0.45),
      q55 = quantile(value, 0.55),
      .groups = "drop"
    ) %>%
    mutate(geography = if_else(
             !is.na(variant),
             as.character(variable),
             NA_character_
           ),
           variable = if_else(
             is.na(variant),
             as.character(variable),
             variant))
  data <- data %>%
    bind_rows(cum_data)
}

data <- data %>%
  mutate(variable = if_else(variable == "2-10", "02-10", variable),
         name = if_else(name == "R", "Rt", name))  

var_names <- c(
  est_prev = "Prevalence estimate",
  infections = "Estimated incidence",
  dcases = "Daily estimated prevalence",
  r = "Growth rate",
  Rt = "Reproduction number",
  cumulative_infections = "Cumulative attack rate"
)

labels <- list(
  est_prev = scales::percent_format(1L, scale = 1),
  infections = scales::comma,
  dcases = scales::percent_format(1L),
  r = waiver(),
  Rt = waiver(),
  cumulative_infections = scales::percent_format(1L)
)

group <- c(
  age_school = "Age group",
  national = "Nation",
  regional = "Region",
  variant_regional = "Variant",
  variant_national = "Variant"
)

non_variant <- data %>%
  filter(!grepl("^variant", level))

histories <- list(all = Inf, `1year` = months(12), `3months` = months(3))
breaks <- c(all = "4 months", `1year` = "3 month", `3months` = "1 month")

for (history in names(histories)) {
  for (level in unique(non_variant$level)) {
    colour_var <- ifelse(level == "local", "region", "variable")
    for (name in unique(non_variant$name)) {
      p <- ggplot(non_variant %>%
                  filter(level == {{ level }},
                         name == {{ name }},
                         date > max(date) - histories[[history]]),
                  aes_string(x = "date",
			     colour = colour_var, fill = colour_var,
			     group = "variable")) +
        geom_ribbon(mapping = aes(ymin = `q45`, ymax = `q55`), alpha = 0.5) +
        geom_ribbon(mapping = aes(ymin = `q25`, ymax = `q75`), alpha = 0.25) +
        geom_ribbon(mapping = aes(ymin = `q5`, ymax = `q95`), alpha = 0.125) +
        ylab(name) + xlab("") +
        scale_x_date(breaks = breaks[[history]],
                     labels = date_format("%b %Y")) +
        scale_y_continuous(var_names[name], labels = labels[[name]]) +
        scale_colour_brewer(group[level],  palette = "Set1") +
        scale_fill_brewer(group[level], palette = "Set1") +
        theme_minimal()
      ggsave(here::here("figures", "additional",
			paste0(level, "_", name, "_",
                               history, ".svg")), p,
             width = 10, height = 5)
    }
  }
}

variant <- data %>%
  filter(grepl("^variant", level))

for (history in names(histories)) {
  for (level in unique(variant$level)) {
    for (name in unique(variant$name)) {
      p <- ggplot(
        variant %>%
          filter(
            level == {{ level }},
            name == {{ name }}
          ),
        aes(x = date, colour = variable, fill = variable)
      ) +
        geom_ribbon(mapping = aes(ymin = `q45`, ymax = `q55`), alpha = 0.5) +
        geom_ribbon(mapping = aes(ymin = `q25`, ymax = `q75`), alpha = 0.25) +
        geom_ribbon(mapping = aes(ymin = `q5`, ymax = `q95`), alpha = 0.125) +
        ylab(name) +
        xlab("") +
        scale_x_date(
          breaks = breaks[[history]],
          labels = date_format("%b %Y")
        ) +
        scale_y_continuous(var_names[name], labels = labels[[name]]) +
        scale_colour_brewer(group[level], palette = "Set1") +
        scale_fill_brewer(group[level], palette = "Set1") +
        theme_minimal() +
        facet_wrap(~geography)
      ggsave(here::here("figures", "additional",
			paste0(level, "_", name, ".svg")), p,
        width = 14, height = 8
      )
    }
  }
}

## separate plots

level <- "national"
name <- "est_prev"
df_plot <- variant %>%
  filter(
    level == {{ level }},
    name == {{ name }}
  ) %>%
  pivot_longer(starts_with("q"), names_to = "quantile") %>%
  mutate(value = value / 100) %>%
  pivot_wider(names_from = "quantile")

p <- ggplot(df_plot %>% filter(variable == "England"),
  aes(x = date, colour = variable, fill = variable)
) +
  geom_ribbon(mapping = aes(ymin = `q45`, ymax = `q55`), alpha = 0.5) +
  geom_ribbon(mapping = aes(ymin = `q25`, ymax = `q75`), alpha = 0.25) +
  geom_ribbon(mapping = aes(ymin = `q5`, ymax = `q95`), alpha = 0.125) +
  ylab(name) +
  xlab("") +
  expand_limits(x = as.Date("2020-03-01")) +
  scale_x_date(breaks = "4 months", labels = date_format("%b %Y")) +
  scale_y_continuous("Percentage testing positive",
                     labels = scales::percent_format(1L)) +
  scale_colour_brewer("Age group", palette = "Set1") +
  scale_fill_brewer("Age group", palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "none")

level <- "age_school"
name <- "est_prev"
df_plot <- variant %>%
  filter(
    level == {{ level }},
    name == {{ name }},
    date >"2021-05-01"
  ) %>%
  pivot_longer(starts_with("q"), names_to = "quantile") %>%
  mutate(value = value / 100) %>%
  pivot_wider(names_from = "quantile")

p <- ggplot(df_plot,
  aes(x = date, colour = variable, fill = variable)
) +
  geom_ribbon(mapping = aes(ymin = `q45`, ymax = `q55`), alpha = 0.5) +
  geom_ribbon(mapping = aes(ymin = `q25`, ymax = `q75`), alpha = 0.25) +
  geom_ribbon(mapping = aes(ymin = `q10`, ymax = `q95`), alpha = 0.125) +
  ylab(name) +
  xlab("") +
  expand_limits(x = as.Date("2021-04-01")) +
  scale_x_date(breaks = "2 months", labels = date_format("%b %Y")) +
  scale_y_continuous("Percentage testing positive",
                     labels = scales::percent_format(1L)) +
  scale_colour_brewer("Age group", palette = "Dark2") +
  scale_fill_brewer("Age group", palette = "Dark2") +
  theme_minimal() +
  ggtitle("England, by age") +
  theme(legend.position = "bottom")
