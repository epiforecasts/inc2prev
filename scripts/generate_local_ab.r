## load libraries
library("here")
library("data.table")
library("tidyr")
library("inc2prev")

## load scripts in `R/`
source(here::here("scripts", "read.R"))

## load vaccination data
vacc <- data.table(read_vacc())[level == "local"]
## load ONS/LTLA/LAD/region mapping
areas <- fread(here::here("data-processed", "cis_areas.csv"))
## England only for now
areas <- areas[!(region %in% c("Wales", "Scotland", "Northern Ireland"))]

## load samples from local prealence model
local_samples <- readRDS(here::here("outputs", "samples_local.rds"))
## load data used for local prevalence model
local_diag <- readRDS(here::here("outputs", "diagnostics_local.rds"))
local_dat <- local_diag$data[[1]]
## load samples from regional antibody model
regional_samples <- readRDS(here::here("outputs", "samples_regional_ab_higher.rds"))
## load data used for local regional antibody model
regional_diag <- readRDS(here::here("outputs", "diagnostics_regional_ab_higher.rds"))
regional_dat <- regional_diag$data[[1]]

## list of parameters to grab from each model
prev_params <- c("alpha", "rho", "eta", "init_inc", "prob_detect", "sigma")
ab_params <- c("beta", "gamma", "delta", "k", "l", "init_dab", "ab_sigma", "ab_sd2")

## list of data sets to grab from each model
dat <- c(local_dat[c("M", "L", "t", "diff_order", "prev_stime", "prev_etime", "ut", "obs")], 
	 regional_dat[c("vacc_ab_delay", "inf_ab_delay")])

## get modelled infection dates for later
inf_dates <- unique(
  local_samples[variable %in% areas$geography_code & name == c("infections")]$date
)

## load in local vaccination data, fill zeroes and convert to matrix
vacc_dates <- unique(regional_samples[name == c("dab")]$date)
vacc_base <- 
  expand.grid(date = vacc_dates,
              variable = na.omit(unique(local_samples$variable)))
vacc <- vacc[, .(
  date = as.Date(date),
  vaccinated = vaccinated,
  variable
)]
setkey(vacc, date, variable)
vacc <- vacc[J(vacc_base), roll = 0]
vacc <- vacc[is.na(vaccinated), vaccinated := 0]
dat[["vacc"]] <- t(as.matrix(
  dcast(vacc, date ~ variable, value.var = "vaccinated")[, -1]
))
dat[["n"]] <- dat[["n_ab"]] <- nrow(dat[["vacc"]])
dat[["ab_index"]] <- seq_len(dat[["n"]])

## function to combine parameters from local prevalence and regional antibody model
combine_params <- function(prev, ab, transfer = "init_dab") {
  for (loop_transfer in transfer) {
    prev[[1]] <- rbind(
      prev[[1]][name != loop_transfer],
      merge(
        areas[, list(variable = geography_code, region, n_index = 1:.N)],
        ab[[1]][name == loop_transfer, list(region = variable, lp__, sigma, name, value, t_index, p_index = NA_integer_, date, level)],
        by = "region"
      )[, -c("region"), with = TRUE]
    )
  }
  return(list(rbind(prev[[1]][name %in% c(prev_params, transfer)], 
		    ab[[1]][name %in% ab_params])))
}

## create parameter list and fix if desired
get_params <- function(dt, fix = list()) {
  dt <- dt[[1]]
  ret <- list()
  for (loop_name in unique(dt$name)) {
    name_dt <- dt[name == loop_name]
    if (all(!is.na(name_dt$variable)) && 
	length(unique(name_dt$variable)) > length(unique(name_dt$n_index))) {
      name_dt <- name_dt[, n_index := as.integer(factor(variable))]
    }
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

## wrapper for i2p_simulate
simulate <- function(data, parameters) {
  sim <- data.table(i2p_simulate(data, parameters))
  setnames(sim, as.character(inf_dates))
  sim <- sim[, local_area := rownames(data$vacc)]
  return(sim)
}

## combine samples from the two models
nls <- nest(local_samples, data = !sample)
nrs <- nest(regional_samples, data = !sample)
samples <- merge(nls, nrs, by = "sample", suffixes = c(".local", ".regional"))
## combine parameters to create a single set of parameters to run model with
samples <- samples[, list(data = combine_params(data.local, data.regional)), by = sample]
## convert parameter samples to lists
samples <- samples[, list(param = get_params(data)), by = sample]
## simulate
sims <- samples[, list(sim = list(simulate(dat, param[[1]]))), by = sample]
## tidy
lsims <- melt(data.table(unnest(sims, cols = c(sim))), id.vars = c("local_area", "sample"), variable.name = "date")
lsims <- lsims[, name := "antibodies"]
lsims <- lsims[, date := as.Date(as.character(date))]
linf <- local_samples[name %in% c("infections", "R"), list(local_area = variable, sample, date, value, name)]

inf_ab <- rbind(lsims, linf)

inf_ab <- inf_ab[local_area %in% unique(areas$geography_code)]
saveRDS(as_tibble(inf_ab), here::here("outputs", "inf_ab_local_higher_samples.rds"))
inf_ab_summary <- inf_ab[, as.list(quantile(.SD, prob = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)), by = c("local_area", "date", "name"), .SDcols = c("value")]
percentages <- grep("%$", colnames(inf_ab_summary), value = TRUE)
qs <- paste0("q", as.numeric(sub("%$", "", percentages)))
setnames(inf_ab_summary, percentages, qs)
saveRDS(as_tibble(inf_ab_summary), here::here("outputs", "inf_ab_local_higher_estimates.rds"))
fwrite(inf_ab_summary, here::here("outputs", "inf_ab_local_higher.csv"))
