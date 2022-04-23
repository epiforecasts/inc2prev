library(data.table)
library(socialmixr)
library(here)

ons_age_limits <- c(2, 11, 16, 25, 35, 50, 70)

vacc <- fread("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=vaccinationsAgeDemographics&format=csv")
vacc_local <- fread("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=vaccinationsAgeDemographics&format=csv")
areas <- fread(here::here("data-processed", "cis_areas.csv"))

vacc[, lower_age_limit := as.integer(sub("[_+].*$", "", age))]

vacc[, vaccinated := rowSums(.SD), .SDcols = grep("^newPeopleVaccinated", names(vacc))]
vacc_local[, vaccinated := rowSums(.SD), .SDcols = grep("^newPeopleVaccinated", names(vacc_local))]

## need to merge twice to dealt with local area remapping
vacc_local <- merge(
  vacc_local, 
  areas[, list(geography_code, areaName = lad, ltla_name, region)], 
  by = "areaName", all.x = TRUE
)
vacc_local[, areaName2 := areaName]
vacc_local <- merge(
  vacc_local, 
  unique(areas[, list(geography_code2 = geography_code, areaName2 = dashboard_name, region2 = region)]),
  by = "areaName2", all.x = TRUE
)
vacc_local[is.na(geography_code), c("geography_code", "region") := list(geography_code2, region2)]

vacc[lower_age_limit == min(lower_age_limit), 
     lower_age_limit := max(ons_age_limits[ons_age_limits <= unique(lower_age_limit)])]
vacc[, lower_age_limit := reduce_agegroups(lower_age_limit, ons_age_limits)]
vacc <- vacc[, list(vaccinated = sum(vaccinated)), 
     by = list(vaccination_date = date, areaName, lower_age_limit)]

vacc_national <- vacc[, list(vaccinated = sum(vaccinated)), by = vaccination_date]
vacc_regional <- vacc[, list(vaccinated = sum(vaccinated)), by = list(vaccination_date, areaName)]
vacc_local <- vacc_local[, list(vaccinated = sum(vaccinated)), by = list(vaccination_date = date, geography_code, region)]
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
  vacc_local[, list(level = "local",
		    vaccination_date,
		    vaccinated,
		    geography = geography_code, region,
		    lower_age_limit = NA_integer_)],
  vacc_age[, list(level = "age_school",
		       vaccination_date,
		       vaccinated,
		       geography = "England",
		       lower_age_limit)]
  ), 
  fill = TRUE
)
setkey(vacc_all, level, geography, lower_age_limit, vaccination_date)
fwrite(vacc_all, here::here("data-processed", "vacc.csv"))
