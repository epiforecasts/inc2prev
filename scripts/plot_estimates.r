library(here)
library(purrr)
library(ggplot2)
library(socialmixr)
library(dplyr)
library(readr)
library(scales)
library(tidyr)

## whether to make spaghetti plots
spaghetti <- FALSE

## Get tools
source(here::here("R", "read.R"))
prev <- read_cis()
local_region <- prev %>%
  filter(level == "local") %>%
  select(level, variable, region) %>%
  distinct()

dir.create(here::here("figures", "additional"),
           showWarnings = FALSE, recursive = TRUE)

var_names <- c(
  est_prev = "Prevalence estimate",
  infections = "Estimated incidence",
  dcases = "Daily estimated prevalence",
  r = "Growth rate",
  Rt = "Reproduction number",
  cumulative_infections = "Cumulative attack rate",
  cumulative_exposure = "Cumulative exposure"
)

labels <- list(
  est_prev = scales::percent_format(1L, scale = 1),
  infections = scales::percent_format(0.1),
  dcases = scales::percent_format(1L),
  r = waiver(),
  Rt = waiver(),
  cumulative_infections = scales::percent_format(1L),
  cumulative_exposure = scales::percent_format(1L)
)

group <- c(
  age_school = "Age group",
  national = "Nation",
  regional = "Region",
  variant_regional = "Variant",
  variant_national = "Variant"
)

hline <- c(
  r = 0,
  Rt = 1
)

histories <- list(all = Inf, `1year` = months(12), `3months` = months(3))
breaks <- c(all = "4 months", `1year` = "3 month", `3months` = "1 month")

file_pattern <-
  paste0(ifelse(spaghetti, "samples", "estimates"), "_[^_]+\\.rds")

files <-
  list.files(here("outputs"),
	     pattern = file_pattern,
	     full.names = TRUE)

for (file in files) {
  data <- readRDS(file) %>%
    separate(variable, c("variable", "geography"), sep = "\\|")
  cum_file <-
    here::here("outputs", paste0("cumulative_", sub("^[^_]+_", "", file)))
  if (file.exists(cum_file)) {
    cum_data <- readRDS(cum_file) %>%
      filter(name %in% names(var_names))
    if (any(grepl("^variant", unique(data$level)))) {
      cum_data <- cum_data %>%
        mutate(geography = as.character(variable),
               variable = variant)
    }
    if (!spaghetti) {
      quantiles <- parse_number(grep("^q", colnames(data), value = TRUE)) / 100
      cum_data <- cum_data %>%
        group_by_at(vars(intersect(colnames(data), colnames(.)))) %>%
        summarise(
          value = quantile(value, quantiles),
          q = paste0("q", quantiles * 100),
          .groups = "drop"
        ) %>%
        pivot_wider(names_from = "q")
    }
    data <- data %>%
      bind_rows(cum_data)
    if (!"geography" %in% colnames(data)) {
      data <- data %>%
        mutate(geography = NA_character_,
               variable = as.character(variable))
    }
  }

  for (level in unique(data$level)) {
    level_data <- data %>%
      filter(level == {{ level }}) %>%
      mutate(variable = if_else(variable == "2-10", "02-10", variable),
             name = if_else(name == "R", "Rt", name)) %>%
      filter(name %in% names(var_names))

    if (level == "local") {
      level_data <- level_data %>%
        left_join(local_region, by = c("level", "variable"))
    }

    colour_var <- ifelse(level == "local", "region", "variable")

    for (history in names(histories)) {
      for (name in unique(level_data$name)) {
        plot_df <- level_data %>%
          filter(name == {{ name }}, 
		 date > max(date) - histories[[history]])
        aes_str <- list(x = "date", colour = colour_var)
        if (spaghetti) {
          plot_df <- plot_df %>%
            mutate(var_sample = interaction(variable, sample))
          aes_str <- c(aes_str, list(y = "value", group = "var_sample"))
        } else {
          aes_str <- c(aes_str, list(fill = colour_var))
        }
        p <- ggplot(plot_df, do.call(aes_string, aes_str)) +
          ylab(name) + xlab("")
        if (spaghetti) {
          p <- p +
            geom_line(alpha = 0.25)
        } else {
          p <- p +
            geom_ribbon(mapping = aes(ymin = `q45`, ymax = `q55`), alpha = 0.5) +
            geom_ribbon(mapping = aes(ymin = `q5`, ymax = `q95`), alpha = 0.125)
        }
        p <- p +
          scale_x_date(breaks = breaks[[history]],
                       labels = date_format("%b %Y")) +
          scale_y_continuous(var_names[name], labels = labels[[name]]) +
          scale_colour_brewer(group[level],  palette = "Set1") +
          scale_fill_brewer(group[level],  palette = "Set1") +
          theme_minimal()
        if (name %in% names(hline)) {
          p <- p +
            geom_hline(yintercept = hline[[name]], linetype = "dashed")
        }
        if (grepl("^variant_", level)) {
          p <- p +
            facet_wrap(~geography)
        }
        ggsave(here::here("figures", "additional",
                          paste0(level, "_", name, "_",
                                 history, ".svg")), p,
               width = 10, height = 5)
      }
    }
  }
}
