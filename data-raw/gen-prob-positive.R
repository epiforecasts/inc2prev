library(boot)
library(data.table)
library(boot)
library(data.table)

dt <- fread("data-raw/hellewell-prI-posterior.csv")

# Can pass a vector of quantiles
prob_detectable <- purrr::map_df(0:60, prob_detectable, dt = dt)
prob_detectable <- data.table::dcast(prob_detectable, sample ~ time,
  value.var = "p"
)

# save
data.table::fwrite(prob_detectable, "data/prob_detectable.csv")
