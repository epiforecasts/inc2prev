#! /usr/bin/env RScript
suppressMessages(library(cmdstanr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(here))
suppressMessages(library(socialmixr))
suppressMessages(library(lubridate))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(future.apply))
suppressMessages(library(future.callr))
suppressMessages(library(future))
suppressMessages(library(docopt))

doc <- "
Estimate incidence from ONS positivity prevalence data,
possibly including antibody and vaccination data
Usage:
    estimate.R

Options:
    -h --help Show this screen
    -a --ab   Use antibody data
"

## if running interactively can set opts to run with options
if (interactive()) {
  if (!exists("opts")) opts <- list()
} else {
  opts <- docopt(doc)
}

antibodies <- !is.null(opts$ab) && opts$ab

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

# Load prevalence data and split by location
data <- read_cis() %>%
  nest(prevalence = c(-variable))

if (antibodies) {
  ab <- read_ab() %>%
    nest(antibodies = c(-variable))
  vacc <- read_vacc() %>%
    nest(vaccination = c(-variable))
  early <- read_early() %>%
    nest(initial_antibodies = c(-variable))
  data <- data %>%
    inner_join(ab, by = "variable") %>%
    inner_join(vacc, by = "variable") %>%
    inner_join(early, by = "variable")
}

data <- data %>%
  group_split(variable)

# Location probability of detection posterior
prob_detect <- read_prob_detectable()

# Compile incidence -> Prevalence model
if (antibodies) {
  mod <- i2p_model("stan/inc2prev_antibodies.stan")
} else {
  mod <- i2p_model()
}

# Compile tune inverse gamma model
tune <- i2p_gp_tune_model()

## Fit model
dir.create(here::here("outputs"), showWarnings = FALSE)

# create a helper function to estimate the model and apply some
# summary statistics
incidence_with_var <- function(data, pb, model, gp_model) {
  message("Fitting model")

  mod <- cmdstanr::cmdstan_model(
    model$stan_file(),
    include_paths = here::here("stan", "functions")
  )
  safe_incidence <- purrr::safely(incidence)

  variables <- c("est_prev", "infections", "dcases", "r", "R")
  prev <- data$prevalence[[1]]
  if (antibodies) {
    variables <- c(variables, "est_ab", "dab", "beta", "gamma", "delta")
    ab <- data$antibodies[[1]]
    vacc <- data$vaccination[[1]]
    init_ab <- data$initial_antibodies[[1]]
  } else {
    ab <- NULL
    vacc <- NULL
    init_ab <- NULL
  }
  fit <- safe_incidence(
    prev = prev,
    ab = ab,
    vacc = vacc,
    init_ab = init_ab,
    variables = variables,
    prob_detect = pb, parallel_chains = 2, iter_warmup = 250,
    chains = 2, model = mod, adapt_delta = 0.9, max_treedepth = 12,
    data_args = list(gp_tune_model = gp_model)
  )

  if (is.null(fit$result)) {
    fit <- data.table::data.table(
      error = list(fit$error)
    )
  } else {
    fit <- fit$result
  }

  level <- unique(data$prevalence[[1]]$level)
  variable <- data$variable

  fit <- fit[, level := level]
  fit <- fit[, variable := variable]

  return(fit)
}

# Run model fits in parallel
plan(callr, workers = future::availableCores())
est <- future_lapply(
  data, incidence_with_var,
  pb = prob_detect,
  model = mod, gp_model = tune, future.seed = TRUE
)

est <- rbindlist(est, use.names = TRUE, fill = TRUE)
# Add summary information to posterior summary and samples
est[, summary := map2(summary, variable, ~ as.data.table(.x)[, variable := .y])]
est[, summary := map2(summary, level, ~ as.data.table(.x)[, level := .y])]
est[, samples := map2(samples, variable, ~ as.data.table(.x)[, variable := .y])]
est[, samples := map2(samples, level, ~ as.data.table(.x)[, level := .y])]

# Bind posterior samples/summary together
estimates <- bind_rows(est$summary)
samples <- bind_rows(est$samples)
diagnostics <- select(est, -samples, -summary)

suffix <- ifelse(antibodies, "_ab", "")

# Save output
saveRDS(samples, paste0("outputs/samples", suffix, ".rds"))
saveRDS(estimates, paste0("outputs/estimates", suffix, ".rds"))
fwrite(estimates, paste0("outputs/estimates", suffix, ".csv"))
saveRDS(diagnostics, paste0("outputs/diagnostics", suffix, ".rds"))
