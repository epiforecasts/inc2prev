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

# Test target
example_var <- "England"
end_date <- "2022-01-01"
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
mod <- i2p_model("stan/inc2prev_antibodies.stan")

# Compile tune inverse gamma model
tune <- i2p_gp_tune_model()

# Fit the infection to prevalence model
fit <- incidence(
  joint_data$prevalence[[1]],
  joint_data$antibodies[[1]],
  joint_data$vaccination[[1]],
  joint_data$initial_antibodies[[1]],
  variables = c(
    "est_prev", "est_ab", "infections", "dcases",
    "dab", "r", "R", "beta", "gamma", "delta"
  ),
  prob_detect = prob_detect, parallel_chains = 2, iter_warmup = 250,
  chains = 2, model = mod, adapt_delta = 0.8, max_treedepth = 12,
  data_args = list(gp_tune_model = tune, ht = 14),
  keep_fit = TRUE
)
fit

# plot modelled and observed (but also modelled) prevalence
prev_plot <- plot_prev(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$prevalence[[1]]
) +
  scale_y_continuous(tran = scales::logit_trans())
ggsave("figures/prev.png", prev_plot, width = 9, height = 6)

# plot modelled and observed (but also modelled) antibodies
ab_plot <- plot_prev(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$antibodies[[1]],
  data_source = "ONS Antibodies", observed = "est_ab",
  modelled = "dab"
)
ggsave("figures/ab.png", ab_plot, width = 9, height = 6)


# pairs plot
stanfit <- read_stan_csv(fit$fit[[1]]$output_files())
np <- nuts_params(stanfit)
pairs <- mcmc_pairs(fit$fit[[1]]$draws(),
  np = np,
  pars = c("alpha", "rho", "beta", "gamma[1]", "gamma[2]", "delta")
)
ggsave("figures/pairs.png", pairs, width = 16, height = 16)

# plot infections
plot_trace(
  fit$samples[[1]][sample <= 100], "infections"
) +
  labs(y = "Infections", x = "Date") +
  scale_y_continuous(labels = scales::percent)
ggsave("figures/infections.png", width = 9, height = 6)

# plot growth
plot_trace(
  fit$samples[[1]][sample <= 100], "r"
) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)
ggsave("figures/growth.png", width = 9, height = 6)

# plot Rt
plot_trace(
  fit$samples[[1]][sample <= 100], "R"
) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)
ggsave("figures/Rt.png", width = 9, height = 6)
