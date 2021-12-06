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
  group_split(variable)

# Location probability of detection posterior
prob_detect <- fread("data/prob_detectable.csv")

# Compile incidence -> Prevalence model
mod <- i2p_model("stan/inc2prev_antibodies.stan")

# Compile tune inverse gamma model
tune <- rstan::stan_model("stan/tune_inv_gamma.stan")

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

  fit <- safe_incidence(
    data$prevalence[[1]],
    data$antibodies[[1]],
    data$vaccination[[1]],
    data$initial_antibodies[[1]],
    variables = c(
      "est_prev", "est_ab", "infections", "dcases",
      "dab", "r", "R", "beta", "gamma", "delta"
    ),
    prob_detect = pb, parallel_chains = 2,
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

  start_date <- min(data$prevalence[[1]]$start_date)
  dates <- data$prevalence[[1]]$date
  fit <-
    fit[, summary := map(
      summary, ~ as.data.table(.x)[
        , date :=
          index2date(
            name, index, start_date,
            dates, data[[1]]$ut
          )
      ]
    )]
  fit <-
    fit[, samples := map(
      samples, ~ as.data.table(.x)[
        , date :=
          index2date(
            name, index, start_date,
            dates, data[[1]]$ut
          )
      ]
    )]

  return(fit)
}

# Run model fits in parallel
plan(callr, workers = future::availableCores())
est <- future_lapply(
  joint_data, incidence_with_var,
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
estimates <- rbindlist(est$summary)
samples <- rbindlist(est$samples)
diagnostics <- est[, -c("samples", "summary"), with = FALSE]

# Save output
saveRDS(samples, "outputs/samples_ab.rds")
saveRDS(estimates, "outputs/estimates_ab.rds")
saveRDS(diagnostics, "outputs/diagnostics_ab.rds")
