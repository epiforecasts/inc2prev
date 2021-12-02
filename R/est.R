## Packages
library(cmdstanr)
library(data.table)
library(EpiNow2)
library(dplyr)
library(ggplot2)
library(here)
library(socialmixr)
library(posterior)
library(readr)

## Get tools
source("R/utils.R")

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
  inner_join(early, by = "variable")

prob_detectable <- fread("data/prob_detectable.csv")

## Format data
dat <- apply(joint_data, 1, function(x) {
  stan_data(x$prev, x$antibodies, x$vaccination, x$initial_antibodies,
            prob_detectable)
})

# Model prep
mod <- cmdstan_model("stan/model_antibodies.stan",
  include_paths = c("stan/functions", "ctdist/stan/functions"),
  cpp_options = list(stan_threads = FALSE)
)
inits <- lapply(dat, stan_inits)

## Fit model
if (!dir.exists(here::here("outputs"))) dir.create(here::here("outputs"))
variables <- c("est_prev", "est_ab", "infections", "dcases", "dab", "r", "R", "beta", "gamma", "delta")
quantiles <- seq(0.05, 0.95, by = 0.05)
nb_samples <- 100

results <- list(estimates = list(), samples = list())
for (i in 1) {
  level <- unique(joint_data[i, ][["prevalence"]][[1]]$level)
  variable <- joint_data[i, ]$variable
  fit <-
    mod$sample(
          data = dat[[i]],
          init = inits[[i]],
          parallel_chains = 2,
          threads_per_chain = 1
        )
  results$estimates[[i]] <-
    fit$summary(variables = variables,
                ~ quantile(.x, probs = quantiles)) %>%
    rename(name = variable) %>%
    mutate(level = level,
           variable = variable)
  results$samples[[i]] <- fit$draws(variables) %>%
    as_draws_df() %>%
    as_tibble() %>%
    mutate(level = level,
           variable = variable) %>%
    mutate(sample = 1:n()) %>%
    filter(sample > n() - nb_samples) %>%
    mutate(sample = 1:n()) %>%
    select(-.chain, -.iteration, -.draw) %>%
    pivot_longer(c(-level, -variable, -sample)) %>%
    pivot_wider(names_from = "sample")
}

## set dates
for (i in 1:length(prev)) {
  for (output in names(results)) {
    start_date <- min(joint_data[i, ][["prevalence"]][[1]]$start_date)
    dates <- joint_data[i, ][["prevalence"]][[1]]$date
    results[[output]][[i]] <- results[[output]][[i]] %>%
      mutate(
        index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
        name = sub("\\[.*$", "", name),
        date = case_when(
          name %in% c("infections", "dcases", "dab") ~
            index - 1 + start_date - dat[[i]]$ut,
          name == "est_prev" ~ dates[index],
          name == "r" ~ index + start_date - dat[[i]]$ut,
          name == "R" ~ index - 1 + start_date
        )
      )
  }
}

saveRDS(bind_rows(results$samples), "outputs/samples.rds")
saveRDS(bind_rows(results$estimates), "outputs/estimates.rds")
