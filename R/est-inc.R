## Packages
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) # nolint
library(cmdstanr)
# install_cmdstan() # nolint
library(data.table)
library(EpiNow2)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(posterior)
library(rstan)
library(truncnorm)
library(purrr)
color_scheme_set("brightblue")

## Get tools
source("R/est-inc-utils.R")

## Read in data
prev <- fread("data/ons-prev.csv")
prob_detectable <- fread("data/prob_detectable.csv")
min_date <- as.Date(min(prev$start_date))

## Format data
region <- "England"
dat <- stan_data(prev, prob_detectable,
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
  pars = c("alpha", "rho", "eta[1]", "eta[2]", "sigma")
)
ggsave("figures/divergent-transitions.png", dts, width = 7, height = 5)

# pairs plot
pairs <- mcmc_pairs(fit$draws(),
  np = np,
  pars = c("alpha", "rho", "eta[1]", "eta[2]", "sigma")
)
pairs
ggsave("figures/pairs.png", pairs, width = 16, height = 16)

## Output

# plot estimated prevalence
plot_prev(fit, prev[geography %in% region])
ggsave("figures/prevalence.png", width = 7, height = 5)

# plot infections
plot_trend(fit, "infections", date_start = min_date - dat$ut) +
  labs(y = "Infections", x = "Date")
ggsave("figures/infections.png", width = 7, height = 5)

# plot population prevalence
plot_trend(fit, "pop_prev", date_start = min_date) +
  labs(y = "Prevalence", x = "Date")
ggsave("figures/population-prevalence.png", width = 7, height = 5)

# plot growth
plot_trend(fit, "r", date_start = min_date - dat$ut - 1) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)
ggsave("figures/growth.png", width = 7, height = 5)

# plot Rt
plot_trend(fit, "R", date_start = min_date - dat$ut + 7) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)
ggsave("figures/Rt.png", width = 7, height = 5)
