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
## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

# Load prevalence data and split by location
prev <- read_cis() %>%
  nest(prevalence = c(-variable))
ab <- read_ab() %>%
  nest(antibodies = c(-variable))
vacc <- read_vacc() %>%
  nest(vaccination = c(-variable))
early <- read_early() %>%
  nest(initial_antibodies = c(-variable))

joint_data <- prev %>%
  inner_join(ab, by = "variable") %>%
  inner_join(vacc, by = "variable") %>%
  inner_join(early, by = "variable") %>%
  filter(variable == example_var) %>%
  group_split(variable)

# Location probability of detection posterior
prob_detect <- fread("data/prob_detectable.csv")

# Compile incidence -> Prevalence model
mod <- i2p_model("stan/inc2prev_antibodies.stan")

# Compile tune inverse gamma model
tune <- rstan::stan_model("stan/tune_inv_gamma.stan")

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
  prob_detect = prob_detect, parallel_chains = 2,
  chains = 2, model = mod, adapt_delta = 0.9, max_treedepth = 12,
  data_args = list(gp_tune_model = tune), keep_fit = TRUE
)
