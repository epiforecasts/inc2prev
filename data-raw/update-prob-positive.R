# Based on: https://doi.org/10.1186/s12916-021-01982-x

library(boot)
library(data.table)
library(boot)
library(data.table)
library(here)

dt <- fread("https://raw.githubusercontent.com/cmmid/pcr-profile/main/fitted_params.csv") # nolint

source(here("R", "prob_detectable.R"))
# Can pass a vector of quantiles
pb <- purrr::map_df(0:40, prob_detectable, dt = dt)
pb <- data.table::dcast(pb, sample ~ time,
  value.var = "p"
)

# save posterior detection traces
data.table::fwrite(pb, "data/prob_detectable_samples.csv")

# summarise probability of detection posterior
pb <- melt(
  copy(pb),
  value.name = "p", id.vars = "sample"
)
pb[, time := as.numeric(as.character(variable))]
pb <- pb[, .(
  median = median(p),
  mean = mean(p),
  sd = sd(p),
  q5 = quantile(p, 0.05),
  q95 = quantile(p, 0.95)
),
  by = time
]
pb <- pb[,
  purrr::map(.SD, signif, digits = 3),
  .SDcols = c("mean", "median", "sd", "q5", "q95"),
  by = time
]

data.table::fwrite(pb, "data/prob_detectable_summary.csv")

# Summarise model posteriors
dt[, p := NULL]
dt <- melt(dt)
dt <- dt[, .(mean = mean(value), sd = sd(value)), by = c("variable")]
dt <- rbind(
  dt,
  data.table(
    variable = c("time", "mean_pb"),
    mean = c(40, mean(pb$mean))
  ),
  fill = TRUE
)

dt <- rbind(
  dt,
  pb[, .(variable = paste0("prob_detect", time), mean, sd)]
)

data.table::fwrite(dt, "data/prob_detectable_params.csv")
