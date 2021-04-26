## Packages
library(cmdstanr)
library(data.table)
library(EpiNow2)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(posterior)
library(rstan)
library(truncnorm)
library(purrr)
library(scales)
cols <- color_scheme_set("brightblue")

## Get tools
source("R/utils.R")

## Read in data
prev <- fread("data/ons-prev.csv")
prob_detectable <- fread("data/prob_detectable.csv")
min_date <- as.Date(min(prev$start_date))

## Format data
region <- "England"
dat <- stan_data(prev,
  copy(prob_detectable)[sample <= 100],
  region = region,
  population = 56286961
)

## Model prep
mod <- cmdstan_model("stan/model.stan",
  include_paths = c("stan/functions", "ctdist/stan/functions"),
  cpp_options = list(stan_threads = FALSE)
)

inits <- stan_inits(dat)

## Fit model
fit <- mod$sample(
  data = dat,
  init = inits,
  parallel_chains = 4,
  threads_per_chain = 1
)

## Fit diagnostics
fit$cmdstan_diagnose()

# get posterior samples
draws <- fit$draws()
draws <- as_draws_df(draws)

# get fit as stanfit object
stanfit <- read_stan_csv(fit$output_files())
np <- nuts_params(stanfit)

# plot dts
dts <- mcmc_parcoord(fit$draws(),
  np = np,
  pars = c("alpha", "rho", "eta[1]", "prob_detect[58]", "sigma")
)
ggsave("figures/divergent-transitions.png", dts, width = 7, height = 5)

# pairs plot
pairs <- mcmc_pairs(fit$draws(),
  np = np,
  pars = c("alpha", "rho", "eta[1]", "prob_detect[58]", "sigma")
)
pairs
ggsave("figures/pairs.png", pairs, width = 16, height = 16)

## Output

# plot estimated prevalence
plot_prev(fit, prev[geography %in% region], date_start = min_date, alpha = 0.03)
ggsave("figures/prevalence.png", width = 9, height = 6)

plot_trace(fit, "prob_detect", date_start = 0, rev_time = TRUE) +
  labs(x = "Days since infection", y = "Probability of detection")
ggsave("figures/probability-detection.png", width = 9, height = 6)

# plot infections
plot_trace(fit, "infections", date_start = min_date - dat$ut) +
  labs(y = "Infections", x = "Date") +
  scale_y_continuous(labels = scales::comma)
ggsave("figures/infections.png", width = 9, height = 6)

# plot growth
plot_trace(fit, "r", date_start = min_date - dat$ut - 1) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)
ggsave("figures/growth.png", width = 9, height = 6)

# plot Rt
plot_trace(fit, "R", date_start = min_date - dat$ut + 7) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)
ggsave("figures/Rt.png", width = 9, height = 6)
