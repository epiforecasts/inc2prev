## Packages
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
library(future.apply)
library(progressr)
library(future)

cols <- color_scheme_set("brightblue")

## Get tools
source("R/utils.R")
source("R/incidence.R")

## Read in data
prev <- fread("data/ons-prev.csv")
load("data/prob_detectable.rda")
min_date <- as.Date(min(prev$start_date))

## Format data
dat_list <- incidence_data(copy(prev)[geography %in% "England"],
  copy(prob_detectable)[sample <= 10],
  population = 56286961
)

## Model prep
mod <- stan_model("stan/inc2prev.stan")
options(mc.cores = 4)
fit <- incidence(dat_list[[1]], mod, iter = 1000, warmup = 500, chains = 4)

# nolint start
# plan("multicore")
# with_progress({
#  fits <- incidence_lapply(dat_list,
#    model = mod, iter = 500,
#    warmup = 250, chains = 2, refresh = -1
#  )
# })
## Fit diagnostics
# fit <- combine_incidence_fits(fits)
# nolint end

draws <- spread_draws(fit)
draws <- as.array(fit)
np <- nuts_params(fit)

# plot dts
dts <- mcmc_parcoord(draws,
  np = np,
  pars = c("alpha", "rho", "eta[1]", "eta[10]", "sigma")
)
ggsave("figures/divergent-transitions.png", dts, width = 7, height = 5)

# pairs plot
pairs <- mcmc_pairs(draws,
  np = np,
  pars = c("alpha", "rho", "eta[1]", "eta[10]", "sigma")
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
