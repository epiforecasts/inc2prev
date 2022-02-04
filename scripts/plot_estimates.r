library(here)
library(purrr)
library(ggplot2)
library(socialmixr)
library(dplyr)
library(readr)
library(scales)

## Get tools
files <-
  list.files(here("outputs"), pattern = "estimates.*\\.csv", full.names = TRUE)
data <- map(files, read_csv) %>%
  bind_rows()

cum_files <-
  list.files(here("outputs"), pattern = "cumulative.*\\.rds", full.names = TRUE)
cum_data <- map(cum_files, readRDS) %>%
  bind_rows() %>%
  group_by(name, date, variable, level) %>%
  summarise(`5%` = quantile(value, 0.05),
            `95%` = quantile(value, 0.95),
            `25%` = quantile(value, 0.25),
            `75%` = quantile(value, 0.75),
            `45%` = quantile(value, 0.45),
            `55%` = quantile(value, 0.55),
            .groups = "drop")

all <- data %>%
  bind_rows(cum_data) %>%
  mutate(variable = if_else(variable == "2-10", "02-10", variable))

for (level in unique(data$level)) {
  for (name in unique(data$name)) {
    p <- ggplot(all %>%
         filter(level == {{ level }},
                name == {{ name }}),
         aes(x = date, colour = variable, fill = variable)) +
      geom_ribbon(mapping = aes(ymin = `45%`, ymax = `55%`), alpha = 0.5) +
      geom_ribbon(mapping = aes(ymin = `25%`, ymax = `75%`), alpha = 0.25) +
      geom_ribbon(mapping = aes(ymin = `5%`, ymax = `95%`), alpha = 0.125) +
      ylab(name) + xlab("") +
      scale_x_date(breaks = "2 months", labels = date_format("%b %Y")) +
      scale_y_continuous("Cumulative infections", labels = scales::percent_format(1L)) +
      scale_colour_brewer("Age group",  palette = "Set1") +
      scale_fill_brewer("Age group", palette = "Set1") +
      theme_minimal()
    ggsave(here::here("figures", paste0(level, "_", name, ".pdf")), p, width = 14, height = 8)
  }
}

