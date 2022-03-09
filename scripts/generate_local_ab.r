library("here")
library("data.table")
library("tidyr")

samples <- readRDS(here::here("outputs", "samples_regional_ab.rds"))
diag <- readRDS(here::here("outputs", "diagnostics_regional_ab.rds"))
dat <- diag$data[[1]]

## create param and fix if desired
get_params <- function(dt, fix = list()) {
  dt <- dt[[1]]
  ret <- list()
  for (loop_name in unique(dt$name)) {
    name_dt <- dt[name == loop_name]
    indices <- 
      vapply(grep("_index$", colnames(name_dt), value = TRUE),
             function(type) {
	       ifelse(any(is.na(name_dt[[type]])), 0, max(name_dt[[type]]))
	     }, 0)
    present_indices <- names(indices)[indices > 0]
    ret[[loop_name]] <- array(name_dt$value, dim = indices[present_indices])
  }

  for (loop_name in names(fix)) {
    ret[[loop_name]][] <- fix[[loop_name]]
  }
    
  return(list(ret))
}

samples <- nest(samples, data = !sample)
samples <- samples[, list(param = get_params(data)), by = sample]
sims <- samples[, list(sim = i2p_simulate(dat, param[[1]])), by = sample]

