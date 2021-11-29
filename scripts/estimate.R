## Packages
library(cmdstanr)
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(here)
library(socialmixr)
library(future.apply)
library(future.callr)
library(future)

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

# Load prevalence data and split by location
prev_list <- read_cis() %>%
  group_split(variable)

# Location probability of detection posterior
prob_detect <- fread("data/prob_detectable.csv")

# Compile incidence -> Prevalence model
mod <- i2p_model()
# Compile tune inverse gamma model
tune <- rstan::stan_model("stan/tune_inv_gamma.stan")

## Fit model
dir.create(here::here("outputs"), showWarnings = FALSE)

# create a helper function to estimate the model and apply some
# summary statistics
incidence_with_var <- function(prev) {
  message("Fitting model for: ", unique(prev$variable))

  fit <- incidence(
    prev,
    prob_detect = prob_detect, parallel_chains = 1,
    chains = 2, model = mod, adapt_delta = 0.95, max_treedepth = 15,
    refresh = 0
  )
  fit[, level := unique(prev$level)]
  fit[, variable := unique(prev$variable)]
  return(fit)
}

# Run model fits in parallel
plan(callr, workers = future::availableCores())
est <- future_lapply(prev_list, incidence_with_var, future.seed = TRUE)
est <- rbindlist(est)
# Add summary information to posterior summary and samples
est[, summary := map2(summary, variable, ~ as.data.table(.x)[, variable := .y])]
est[, summary := map2(summary, level, ~ as.data.table(.x)[, level := .y])]
est[, samples := map2(samples, variable, ~ as.data.table(.x)[, variable := .y])]
est[, samples := map2(samples, level, ~ as.data.table(.x)[, level := .y])]

# Bind posterior samples/summary together
estimates <- bind_rows(est$summary)
samples <- bind_rows(est$samples)
diagnostics <- select(est, -samples, -summary)

# Save output
saveRDS(samples, "outputs/samples.rds")
saveRDS(estimates, "outputs/estimates.rds")
saveRDS(diagnostics, "outputs/diagnostics.rds")
