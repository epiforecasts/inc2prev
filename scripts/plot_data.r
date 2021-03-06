library(here)
library(purrr)
library(ggplot2)
library(socialmixr)
library(dplyr)

source(here::here("scripts", "read.R"))
prev <- read_cis()

p <- ggplot(prev %>%
         filter(level == "national"),
       aes(x = end_date, ymin = lower, ymax = upper)) +
  geom_ribbon() +
  scale_x_date(breaks = "2 months", labels = date_format("%b %Y")) +
  theme_minimal() +
  facet_wrap(~ variable) +
  xlab("")

p_natural <- p +
  scale_y_continuous("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

p_log <- p +
  scale_y_log10("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

ggsave(here::here("figures", "national_prevalence.pdf"), p_natural,
       width = 12, height = 8)
ggsave(here::here("figures", "national_prevalence_log.pdf"), p_log,
       width = 12, height = 8)


p <- ggplot(prev %>%
         filter(level == "regional"),
       aes(x = end_date, ymin = lower, ymax = upper)) +
  geom_ribbon() +
  scale_x_date(breaks = "2 months", labels = date_format("%b %Y")) +
  theme_minimal() +
  facet_wrap(~ variable) +
  xlab("")

p_natural <- p +
  scale_y_continuous("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

p_log <- p +
  scale_y_log10("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

ggsave(here::here("figures", "regional_prevalence.pdf"), p_natural,
       width = 16, height = 12)
ggsave(here::here("figures", "regional_prevalence_log.pdf"), p_log,
       width = 16, height = 12)


age <- prev %>%
  filter(level == "age_school") %>%
  mutate(variable =
           factor(variable, levels = c("2-10", "11-15", "16-24", "25-34",
                                       "35-49", "50-69", "70+")))

p <- ggplot(age,
            aes(x = end_date, ymin = lower, ymax = upper,
                colour = variable, fill = variable)) +
  geom_ribbon(alpha = 0.5) +
  scale_x_date(breaks = "2 months", labels = date_format("%b %Y")) +
  scale_colour_brewer("Age", palette = "Set1") +
  scale_fill_brewer("Age", palette = "Set1") +
  theme_minimal() +
  ggtitle("England") +
  xlab("")

p_natural <- p +
  scale_y_continuous("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

p_log <- p +
  scale_y_log10("Prevalence of SARS-CoV-2 positivity",
                label = scales::percent_format(accuracy = 0.1))

ggsave(here::here("figures", "age_prevalence.pdf"), p_natural,
       width = 12, height = 8)
ggsave(here::here("figures", "age_prevalence_log.pdf"), p_log,
       width = 12, height = 8)
