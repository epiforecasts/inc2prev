# packages
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) # nolint
library(cmdstanr)
# install_cmdstan() # nolint
library(data.table)
library(EpiNow2)
library(dplyr)
library(ggplot2)
library(bayesplot)
color_scheme_set("brightblue")


# get tools
source("R/est-inc-utils.R")

# read in data
prev <- fread("data/ons-prev.csv")
prob_detectable <- fread("data/prob_detectable.csv")
min_date <- min(prev$date)

# format data for fitting
dat <- stan_data(prev, prob_detectable)

inits <- stan_inits(dat)

# load model
mod <- cmdstan_model("stan/model.stan",
  include_paths = c("stan/functions", "ctdist/stan/functions"),
  cpp_options = list(stan_threads = TRUE)
)

# fit model
fit <- mod$sample(
  data = dat,
  init = inits,
  parallel_chains = 4,
  threads_per_chain = 1
)

# check
fit$cmdstan_diagnose()

# summarise fit
fit$cmdstan_summary()

# plot infections
plot_trend(fit, "infections", date_start = min_date - dat$ut) +
  labs(y = "Infections", x = "Date")

ggsave("figures/infections.png", width = 7, height = 5)

plot_trend(fit, "r", date_start = min_date - dat$ut - 1) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)

ggsave("figures/growth.png", width = 7, height = 5)

plot_trend(fit, "R", date_start = min_date - dat$ut + 7) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)

ggsave("figures/Rt.png", width = 7, height = 5)
