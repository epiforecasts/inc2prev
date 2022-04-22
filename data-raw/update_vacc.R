library(data.table)
library(socialmixr)
library(here)

ons_age_limits <- c(2, 11, 16, 25, 35, 50, 70)

vacc <- fread("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=vaccinationsAgeDemographics&format=csv")

vacc[, lower_age_limit := as.integer(sub("[_+].*$", "", age))]

vacc[, vaccinated := rowSums(.SD), .SDcols = grep("^newPeopleVaccinated", names(vacc))]
vacc[lower_age_limit == min(lower_age_limit), 
     lower_age_limit := max(ons_age_limits[ons_age_limits <= unique(lower_age_limit)])]
vacc[, lower_age_limit := reduce_agegroups(lower_age_limit, ons_age_limits)]
vacc <- vacc[, list(vaccinated = sum(vaccinated)), 
     by = list(vaccination_date = date, areaName, lower_age_limit)]

vacc_national <- vacc[, list(vaccinated = sum(vaccinated)), by = vaccination_date]
vacc_regional <- vacc[, list(vaccinated = sum(vaccinated)), by = list(vaccination_date, areaName)]
vacc_age <- vacc[, list(vaccinated = sum(vaccinated)), by = list(vaccination_date, lower_age_limit)]

vacc_all <- rbindlist(list(
  vacc_national[, list(level = "national",
		       vaccination_date,
		       vaccinated,
		       geography = "England",
		       lower_age_limit = NA_integer_)],
  vacc_regional[, list(level = "regional",
		       vaccination_date,
		       vaccinated,
		       geography = areaName,
		       lower_age_limit = NA_integer_)],
  vacc_age[, list(level = "age_school",
		       vaccination_date,
		       vaccinated,
		       geography = "England",
		       lower_age_limit)]
))
setkey(vacc_all, level, geography, lower_age_limit, vaccination_date)
fwrite(vacc_all, here::here("data-processed", "vacc.csv"))
