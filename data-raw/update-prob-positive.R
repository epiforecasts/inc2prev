library(boot)
library(data.table)
library(boot)
library(data.table)
library(here)

dt <- fread("https://raw.githubusercontent.com/cmmid/pcr-profile/main/fitted_params.csv") # nolint

source(here("R", "prob_detectable.R"))
# Can pass a vector of quantiles
pb <- purrr::map_df(0:60, prob_detectable, dt = dt)
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
  sd = sd(p)
),
  by = time
]
pb <- pb[,
  purrr::map(.SD, signif, digits = 3),
  .SDcols = c("mean", "median", "sd"),
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
    mean = c(60, mean(pb$mean))
  ),
  fill = TRUE
)

data.table::fwrite(dt, "data/prob_detectable_params.csv")
