# packages
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) # nolint
library(cmdstanr)
# install_cmdstan() # nolint
library(data.table)
library(EpiNow2)

# get tools
source("R/est-inc-utils.R")

# read in data
prev <- fread("data/ons-prev.csv")
prob_detectable <- fread("data/prob_detectable.csv")

# format data for fitting
dat <- stan_data(prev, prob_detectable)

# load model
mod <- cmdstan_model("stan/model.stan",
  include_paths = c("stan", "ctdist/stan/functions"),
  cpp_options = list(stan_threads = TRUE)
)

# fit model
fit <- mod$sample(
  data = stan_data, parallel_chains = 4,
  threads_per_chain = 1
)

# check
fit$cmdstan_diagnose()

# summarise fit
fit$cmdstan_summary()
