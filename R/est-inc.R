# packages
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) # nolint
library(cmdstanr)
# install_cmdstan() # nolint
library(data.table)

# read in data
prev <- fread("data/ons-prev.csv")
prob_detectable <- fread("data/prob_detectable.csv")

# extract a single region for prevalence and build features
region <- "England"
prev <- prev[geography %in% region][, .(date, prev = middle)]
prev[, time := date - min(date)]

# summarise prob_detectable for simplicity
prob_detectable <- melt(
  prob_detectable,
  value.name = "p", id.vars = "sample"
)
prob_detectable <- prob_detectable[, .(p = median(p)), by = variable]
prob_detectable[, time := as.numeric(as.character(variable))]
# initial unobserved time
ut <- 14

# build stan data
stan_data <- list(
  ut = ut,
  ot = max(prev$time),
  t = ut + max(prev$time),
  prev = prev$prev,
  prev_time = prev$time,
  prob_detect = prob_detectable$p,
  prob_detect_time = prob_detectable$time
)
