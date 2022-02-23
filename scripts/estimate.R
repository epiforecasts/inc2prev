#! /usr/bin/env Rscript

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
    estimate.R [--ab] [--local | --regional | --age | --variants] [--nhse] [--differencing=<level>] [--start-date=<date>] [--gp-frac=<frac>]
    estimate.R -h | --help

Options:
    -h, --help                 Show this screen
    -a, --ab                   Use antibody data
    -r, --regional             Model regional dynamics
    -l, --local                Model local dynamics
    -g, --age                  Model age
    -v, --variants             Model variants
    -n, --nhse                 Analyse NHSE regions
    -d, --differencing=<level> Level of differencing of GP (0 = infections,  1 = growth, 2 = differences in growth etc.)
    -e, --start_date=<date>    Start date to use for estimation
    -w, --weekly               Aggregate all data to weekly
    -c, --gp_frac=<frac>       Fraction of latent timepoints to use in the Gaussian
                               process approximation. Reducing this improves runtimes at
                               the cost of reducing the accuracy of the Gaussian process
                               approximation.
"

## if running interactively can set opts to run with options
if (interactive()) {
  if (!exists("opts")) opts <- list()
} else {
  opts <- docopt(doc)
}

antibodies <- !is.null(opts$ab) && opts$ab
regional <- !is.null(opts$regional) && opts$regional
local <- !is.null(opts$local) && opts$local
age <- !is.null(opts$age) && opts$age
variants <- !is.null(opts$variants) && opts$variants
nhse <- !is.null(opts$nhse) && opts$nhse
differencing <- ifelse(is.null(opts$differencing), 0L, as.integer(opts$differencing))
start_date <- as.Date(opts$start_date)
weekly <- !is.null(opts$weekly) && opts$weekly
gp_frac <- ifelse(is.null(opts$gp_frac), 0.3, as.numeric(opts$gp_frac))

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

# Load prevalence data and split by location
data <- read_cis(nhse_regions = nhse)

if (regional) {
  filter_level <- "regional"
  suffix <- "_regional"
} else if (local) {
  filter_level <- "local"
  suffix <- "_local"
} else if (age) {
  filter_level <- "age_school"
  suffix <- "_age"
} else if (variants) {
  filter_level <- c("variant_national", "variant_regional")
  suffix <- "_variants"
} else {
  filter_level <- "national"
  suffix <- "_national"
}

data <- data %>%
  filter(level %in% filter_level)

filter_opt <- function(data, start_date) {
  if (length(start_date) > 0) {
    data <- data %>%
      filter(end_date >= start_date)
  }
  return(data)
}

data <- data %>%
  filter_opt(start_date) %>%
  nest(prevalence = c(-variable))

if (antibodies) {
  ab <- read_ab(nhse_regions = nhse) %>%
    filter_opt(start_date) %>%
    nest(antibodies = c(-variable))
  vacc <- read_vacc(nhse_regions = nhse) %>%
    filter_opt(start_date) %>%
    nest(vaccination = c(-variable))
  early <- read_early(nhse_regions = nhse) %>%
    filter_opt(start_date) %>%
    nest(initial_antibodies = c(-variable))
  data <- data %>%
    left_join(ab, by = "variable") %>%
    left_join(vacc, by = "variable") %>%
    left_join(early, by = "variable")
}

## determine estimation strategy
if (antibodies) {
  if (regional) {
    ## fit all at once
    data <- data %>%
      mutate(grouping = "all")
  } else if (local) {
    stop("No local antibody estimates are available")
  } else if (age) {
    ## fit age groups without data alongside older ones
    data <- data %>%
      mutate(grouping = variable,
	     grouping = recode(grouping,
			       `2-10` = "16-24",
			       `11-15` = "16-24"))
  } else if (variants) {
    ## doesn't really make sense with antibodies
    stop("Fitting variant antibodies does not really make sense")
  } else {
    ## national: fit separately as time series have different lengths
    data <- data %>%
      mutate(grouping = variable)
  }
} else {
  ## fit separately
  data <- data %>%
    mutate(grouping = variable)
}

if (local && !weekly) {
  warning("Converting everything to weekly as needed for local estimates")
  weekly <- TRUE
}

data <- data %>%
  group_split(grouping, .keep = FALSE)

# Location probability of detection posterior
prob_detect <- read_prob_detectable()

# Compile incidence -> Prevalence model
mod <- i2p_model()

# Compile tune inverse gamma model
tune <- i2p_gp_tune_model()

## Fit model
dir.create(here::here("outputs"), showWarnings = FALSE)

convert_to_weekly <- function(x, ref_end_date, cols = c("middle", "lower", "upper"), aggregate = mean) {
  if (!("end_date" %in% colnames(x))) {
    x <- x %>%
      mutate(start_date = date,
	     end_date = date)
  }
  x %>%
    pivot_longer(all_of(cols)) %>%
    mutate(end_date = end_date + as.integer(ref_end_date - end_date) %% 7) %>%
    group_by(across(c(-value, -start_date, -date))) %>%
    summarise(start_date = min(start_date), value = aggregate(value), .groups = "drop") %>%
    mutate(date = start_date + (end_date - start_date) / 2) %>%
    pivot_wider()
}

# create a helper function to estimate the model and apply some
# summary statistics
incidence_with_var <- function(data, pb, model, gp_model, differencing = 0, weekly = FALSE, ...) {
  message("Fitting model")

  safe_incidence <- purrr::safely(incidence)

  prev <- data %>%
    unnest(prevalence)
  if (weekly) {
    prev <- prev %>%
      convert_to_weekly(max(.$end_date))
  }
  if (antibodies) {
    ab <- data %>%
      unnest(antibodies)
    vacc <- data %>%
      unnest(vaccination)
    init_ab <- data %>%
      unnest(initial_antibodies)
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
    var_col = "variable",
    prob_detect = pb, parallel_chains = 2, iter_warmup = 250,
    chains = 2, model = model, adapt_delta = 0.9, max_treedepth = 12,
    data_args = list(
      gp_tune_model = gp_model, gp_m = gp_frac, differencing = differencing
    ),
    ...
  )

  if (is.null(fit$result)) {
    fit <- data.table::data.table(
      error = list(fit$error)
    )
  } else {
    fit <- fit$result
  }

  level <- unique(data$prevalence[[1]]$level)
  fit <- fit[, level := level]

  return(fit)
}

# Run model fits in parallel
plan(callr, workers = future::availableCores() %/% 2)
est <- future_lapply(
  data, incidence_with_var,
  pb = prob_detect,
  model = mod, gp_model = tune,
  differencing = differencing,
  weekly = weekly,
  future.seed = TRUE
)

est <- rbindlist(est, use.names = TRUE, fill = TRUE)
# Add summary information to posterior summary and samples
est[, summary := map2(summary, level, ~ as.data.table(.x)[, level := .y])]
est[, samples := map2(samples, level, ~ as.data.table(.x)[, level := .y])]

# Bind posterior samples/summary together
estimates <- bind_rows(est$summary)
samples <- bind_rows(est$samples)
est[, variable := data$variable]
diagnostics <- select(est, -samples, -summary)

suffix <- paste0(suffix, ifelse(antibodies, "_ab", ""))

# Save output
saveRDS(samples, paste0("outputs/samples", suffix, ".rds"))
saveRDS(estimates, paste0("outputs/estimates", suffix, ".rds"))
saveRDS(diagnostics, paste0("outputs/diagnostics", suffix, ".rds"))

pop <- data %>%
  bind_rows() %>%
  unnest(prevalence) %>%
  group_by(variable, level) %>%
  summarise(population = unique(population), .groups = "drop")
format_estimates <- estimates %>%
  left_join(pop, by = c("level", "variable")) %>% 
  filter(!is.na(population)) %>%
  pivot_longer(matches("^q[0-9]+$"), names_to = "quantile") %>%
  mutate(value = if_else(name == "est_prev", value * 100, value),
	 value = if_else(name == "infections", round(value * population), value)) %>%
  pivot_wider(names_from = "quantile")

fwrite(format_estimates, paste0("outputs/estimates", suffix, ".csv"))
