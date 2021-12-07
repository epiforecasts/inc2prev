## Packages
library(cmdstanr)
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(here)
library(socialmixr)
library(lubridate)
library(readr)
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
  prob_detect = prob_detect, parallel_chains = 2, iter_warmup = 500,
  chains = 2, model = mod, adapt_delta = 0.8, max_treedepth = 12,
  data_args = list(gp_tune_model = tune), keep_fit = TRUE
)
fit

# plot modelled and observed (but also modelled) prevalence
plot_prev(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$prevalence[[1]]
)

# plot modelled and observed (but also modelled) antibodies
plot_ab(
  fit$summary[[1]], fit$samples[[1]][sample <= 100],
  joint_data$antibodies[[1]]
)
