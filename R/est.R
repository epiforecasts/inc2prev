## Packages
library(cmdstanr)
library(data.table)
library(EpiNow2)
library(dplyr)
library(ggplot2)
library(here)
library(socialmixr)

## Get tools
source("R/utils.R")

prev <- read_cis() %>%
  group_split(variable)

prob_detectable <- fread("data/prob_detectable.csv")

## Format data
dat <- lapply(prev, stan_data, prob_detectable)

# Model prep
mod <- cmdstan_model("stan/model.stan",
  include_paths = c("stan/functions"),
  cpp_options = list(stan_threads = FALSE)
)

inits <- lapply(dat, stan_inits)

## Fit model
if (!dir.exists(here::here("outputs"))) dir.create(here::here("outputs"))
variables <- c("pop_prev", "est_prev", "infections", "cumulative_infections", "r", "R")
quantiles <- seq(0.05, 0.95, by = 0.05)
nb_samples <- 100

estimates <- list()
samples <- list()
for (i in 1:length(prev)) {
  fit <-
    mod$sample(
          data = dat[[i]],
          init = inits[[i]],
          parallel_chains = 4,
          threads_per_chain = 1
        )
  estimates[[i]] <- fit$summary(variables = variables, ~ quantile(.x, probs = quantiles)) %>%
	  rename(name = variable) %>%
	  mutate(level = unique(prev[[i]]$level),
		 variable = unique(prev[[i]]$variable))
  samples[[i]] <- fit$draws(variables) %>%
	  as_draws_df() %>%
	  as_tibble() %>%
	  mutate(level = unique(prev[[i]]$level),
		 variable = unique(prev[[i]]$variable)) %>%
	  mutate(sample = 1:n()) %>%
	  filter(sample <= nb_samples)
}

estimates <- bind_rows(estimates) %>%
  mutate(time = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
	       name = sub("\\[.*$", "", name))

samples <- bind_rows(samples) %>%
	select(-.chain, -.iteration, -.draw) %>%
        pivot_longer(matches("[0-9]")) %>%
	mutate(time = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
	       name = sub("\\[.*$", "", name)) %>%
	pivot_wider()

saveRDS(samples, "outputs/samples.rds")
saveRDS(estimates, "outputs/estimates.rds")
