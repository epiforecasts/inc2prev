## Packages
library(bayesplot)
library(cmdstanr)
library(data.table)
library(dplyr)
library(purrr)
library(posterior)
library(ggplot2)
library(here)
library(socialmixr)
library(lubridate)
library(readr)
library(rstan)
library(tidyr)
library(future.apply)
library(future.callr)
library(future)
library(cowplot)

# Test target
example_var <- "England"
end_date <- "2021-11-01"
## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

# Load prevalence data and split by location
prev <- read_cis() %>%
  filter(end_date < {{ end_date }}) %>%
  nest(prevalence = c(-variable))
ab <- read_ab() %>%
  filter(end_date < {{ end_date }}) %>%
  nest(antibodies = c(-variable))
vacc <- read_vacc() %>%
  filter(date < {{ end_date }}) %>%
  nest(vaccination = c(-variable))
early <- read_early() %>%
  nest(initial_antibodies = c(-variable))

joint_data <- prev %>%
  inner_join(ab, by = "variable") %>%
  inner_join(vacc, by = "variable") %>%
  inner_join(early, by = "variable") %>%
  filter(variable == example_var)

# Location probability of detection posterior
prob_detect <- read_prob_detectable()

# Compile incidence -> Prevalence model
mod <- i2p_model()

# Compile tune inverse gamma model
tune <- i2p_gp_tune_model()

# Fit the infection to prevalence model
fit <- incidence(
  joint_data$prevalence[[1]],
  joint_data$antibodies[[1]],
  joint_data$vaccination[[1]],
  joint_data$initial_antibodies[[1]],
  variables = c(
    "est_prev", "infections", "dcases", "r", "R",
    "est_ab", "dab", "beta", "gamma", "delta", "k", "l"
  ),
  prob_detect = prob_detect, parallel_chains = 2, iter_warmup = 200,
  chains = 2, model = mod, adapt_delta = 0.85, max_treedepth = 15,
  data_args = list(
    gp_tune_model = tune, horizon = 14, differencing = 1,
    gp_m = 0.3
  ),
  keep_fit = TRUE
)
dir.create(here::here("outputs"), showWarnings = FALSE)
fit$fit[[1]]$save_object(here::here("outputs", "example-fit.rds"))

dir.create(here::here("figures", "example"), 
	   showWarnings = FALSE, recursive = TRUE)
# plot modelled and observed (but also modelled) prevalence
prev_plot <- plot_prev(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$prevalence[[1]]
) +
  scale_x_date(date_breaks = "4 months", date_label = "%b %Y")
ggsave(here::here("figures", "example", "prev.png"), prev_plot, width = 9, height = 6)

# plot modelled and observed (but also modelled) antibodies
ab_plot <- plot_prev(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$antibodies[[1]],
  data_source = "ONS Antibodies", observed = "est_ab",
  modelled = "dab"
) +
  scale_y_continuous("Antibody prevalence", labels = scales::percent) +
  scale_x_date("Date", date_breaks = "4 months", date_label = "%b %Y")

ggsave(here::here("figures", "example", "ab.png"), ab_plot, width = 9, height = 6)

# pairs plot
stanfit <- read_stan_csv(fit$fit[[1]]$output_files())
np <- nuts_params(stanfit)
pairs <- mcmc_pairs(fit$fit[[1]]$draws(),
  np = np,
  pars = c("beta[1]", "gamma[1]", "gamma[2]", "delta[1]", "k[1]", "l[1]")
)
ggsave(here::here("figures", "example", "pairs.png"), pairs, width = 16, height = 16)

# plot infections
inc_plot <- plot_trace(
  fit$samples[[1]][sample <= 100], "infections"
) +
  scale_y_continuous("Incident infections", labels = scales::percent) +
  scale_x_date("Date", date_breaks = "4 months", date_label = "%b %Y")

ggsave(here::here("figures", "example", "infections.png"), inc_plot, width = 9, height = 6)

# plot growth
growth_plot <- plot_trace(
  fit$samples[[1]][sample <= 100], "r"
) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)
ggsave(here::here("figures", "example", "growth.png"), growth_plot, width = 9, height = 6)

# plot Rt
rt_plot <- plot_trace(
  fit$samples[[1]][sample <= 100], "R"
) +
  labs(y = "Reproduction number") +
  scale_x_date("Date", date_breaks = "4 months", date_label = "%b %Y")
  geom_hline(yintercept = 1, linetype = 2)
ggsave(here::here("figures", "example", "Rt.png"), rt_plot, width = 9, height = 6)

p <- plot_grid(prev_plot, inc_plot, ab_plot, rt_plot, labels = c("A", "B", "C", "D"))
ggsave(here::here("figures", "example", "example-estimates.png"), p, width = 12, height = 6)

params <- fit$summary[[1]][is.na(date)]
params[, c("n_index", "t_index", "date") := NULL]
fwrite(params, "outputs/example-parameters.csv")
