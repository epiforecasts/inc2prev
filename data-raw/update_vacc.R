library(data.table)
library(jsonlite)
library(here)

ons_age_limits <- c(2, 11, 16, 25, 35, 50, 70)

## the coronavirus dashboard API has been retired; its data are archived by UKHSA
archive_url <- paste0(
  "https://archive.ukhsa-dashboard.data.gov.uk/",
  "coronavirus-dashboard/vaccinations.zip"
)
## seasonal booster campaigns are only available from the current dashboard API
api_url <- paste0(
  "https://api.ukhsa-dashboard.data.gov.uk/themes/infectious_disease/",
  "sub_themes/respiratory/topics/COVID-19/geography_types"
)
booster_campaigns <- c("autumn22", "spring23", "autumn23")

raw_dir <- here::here("data-raw", "vacc")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
zip_file <- file.path(raw_dir, "vaccinations.zip")
if (!file.exists(zip_file)) {
  options(timeout = max(3600, getOption("timeout")))
  download.file(archive_url, zip_file, mode = "wb")
}

read_archive <- function(file, ...) {
  fread(
    cmd = paste("unzip -p", shQuote(zip_file), shQuote(file.path("Vaccinations", file))),
    ...
  )
}
archive_files <- unzip(zip_file, list = TRUE)$Name

## first and second doses, and third primary or first booster doses
doses <- c("FirstDose", "SecondDose", "ThirdInjection")

read_doses <- function(area_type) {
  rbindlist(lapply(doses, function(dose) {
    read_archive(
      sprintf("%s_newPeopleVaccinated%sByVaccinationDate.csv", area_type, dose),
      select = c("date", "area_code", "area_name", "value")
    )[, dose := dose]
  }))
}

## national and regional doses in England; local doses cover all survey areas
national <- read_doses("nation")[area_name == "England"]
regional <- read_doses("region")[grepl("^E", area_code)]
local <- read_doses("ltla")

## age distribution of doses in England, from upper- or lower-tier local
## authority breakdowns; neither is complete on every day, so on each day we use
## the more complete one and interpolate the age distribution where both are
## missing, then scale to national totals
age_bands <- c(
  "05_11", "12_15", "16_17", "18_24", "25_29", "30_34", "35_39", "40_44",
  "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75_79", "80_84",
  "85_89", "90+"
)
age_cols <- c(
  FirstDose = "newPeopleVaccinatedFirstDoseByVaccinationDate",
  SecondDose = "newPeopleVaccinatedSecondDoseByVaccinationDate",
  ThirdInjection = "newPeopleVaccinatedThirdDoseByVaccinationDate"
)

read_age <- function(area_type) {
  files <- grep(
    sprintf("^Vaccinations/%s_vaccinationsAgeDemographics", area_type),
    archive_files, value = TRUE
  )
  age <- rbindlist(lapply(basename(files), read_archive,
    select = c("date", "area_code", "ageCategory", unname(age_cols))
  ))
  age <- unique(age[grepl("^E", area_code) & ageCategory %in% age_bands])
  age <- melt(age,
    id.vars = c("date", "area_code", "ageCategory"),
    variable.name = "dose", value.name = "vaccinated"
  )
  age[, dose := names(age_cols)[match(dose, age_cols)]]
  age[, list(
    vaccinated = sum(vaccinated, na.rm = TRUE),
    areas = uniqueN(area_code)
  ), by = list(date, dose, age = ageCategory)]
}

age <- rbind(
  read_age("utla")[, source := "utla"],
  read_age("ltla")[, source := "ltla"]
)
age[, coverage := areas / max(areas), by = source]
## prefer complete upper-tier data, then near-complete lower-tier data
age[, priority := fifelse(source == "utla", 1L, 2L)]
age <- age[coverage >= 0.8]
age <- age[, .SD[priority == min(priority)], by = list(date, dose)]
age[, share := vaccinated / sum(vaccinated), by = list(date, dose)]

all_age_dates <- CJ(
  date = seq(min(national$date), max(national$date), by = "day"),
  dose = doses,
  age = age_bands
)
age <- merge(all_age_dates, age[, list(date, dose, age, share)],
  by = c("date", "dose", "age"), all.x = TRUE
)
age[, share := approx(
  as.numeric(date), share, xout = as.numeric(date), rule = 2
)$y, by = list(dose, age)]
age[is.na(share), share := 0]
age <- merge(age, national[, list(date, dose, total = value)],
  by = c("date", "dose")
)
age[, vaccinated := total * share / sum(share), by = list(date, dose)]
age[, lower_age_limit := as.integer(sub("[_+].*$", "", age))]

## seasonal boosters by region and age
api_get <- function(url) {
  results <- list()
  while (!is.null(url)) {
    page <- fromJSON(url)
    results <- c(results, list(as.data.table(page$results)))
    url <- page[["next"]]
  }
  rbindlist(results)
}

boosters <- rbindlist(lapply(booster_campaigns, function(campaign) {
  metric <- sprintf("COVID-19_vaccinations_%s_dosesByDay", campaign)
  geographies <- list(
    "Nation" = "England",
    "UKHSA Region" = unique(regional$area_name)
  )
  rbindlist(lapply(names(geographies), function(type) {
    rbindlist(lapply(geographies[[type]], function(ons_name) {
      ## UKHSA names differ from ONS region names in one case
      api_name <- sub("The Humber", "Humber", ons_name)
      url <- URLencode(sprintf(
        "%s/%s/geographies/%s/metrics/%s?sex=all&page_size=365",
        api_url, type, api_name, metric
      ))
      api_get(url)[, list(
        date = as.IDate(date), age, vaccinated = metric_value,
        geography_type = type, area_name = ons_name
      )]
    }))
  }))
}))
## keep non-overlapping age bands
boosters <- boosters[grepl("^[0-9]+-[0-9]+$", age) | age == "80+"]
boosters[, lower_age_limit := as.integer(sub("[-+].*$", "", age))]

## combine
map_age_limits <- function(limits) {
  ons_age_limits[pmax(findInterval(limits, ons_age_limits), 1L)]
}

vacc_national <- rbind(
  national[, list(vaccination_date = date, vaccinated = value)],
  boosters[geography_type == "Nation", list(vaccination_date = date, vaccinated)]
)[, list(vaccinated = sum(vaccinated)), by = vaccination_date]

vacc_regional <- rbind(
  regional[, list(vaccination_date = date, geography = area_name, vaccinated = value)],
  boosters[geography_type == "UKHSA Region",
    list(vaccination_date = date, geography = area_name, vaccinated)
  ]
)[, list(vaccinated = sum(vaccinated)), by = list(vaccination_date, geography)]

vacc_age <- rbind(
  age[, list(vaccination_date = date, lower_age_limit, vaccinated)],
  boosters[geography_type == "Nation",
    list(vaccination_date = date, lower_age_limit, vaccinated)
  ]
)[, lower_age_limit := map_age_limits(lower_age_limit)][,
  list(vaccinated = sum(vaccinated)), by = list(vaccination_date, lower_age_limit)
]

## local authorities mapped to infection survey areas; names are matched twice
## to deal with local area remapping, and normalised because the archive
## replaces spaces with hyphens
normalise_name <- function(x) {
  trimws(gsub("[^a-z]+", " ", tolower(x)))
}
areas <- fread(here::here("data-processed", "cis_areas.csv"))
local[, area_name := normalise_name(area_name)]
vacc_local <- merge(
  local, areas[, list(geography_code, area_name = normalise_name(lad), region)],
  by = "area_name", all.x = TRUE
)
vacc_local <- merge(
  vacc_local,
  unique(areas[, list(
    geography_code2 = geography_code,
    area_name = normalise_name(dashboard_name), region2 = region
  )]),
  by = "area_name", all.x = TRUE
)
vacc_local[is.na(geography_code), c("geography_code", "region") :=
  list(geography_code2, region2)]
vacc_local <- vacc_local[!is.na(geography_code),
  list(vaccinated = sum(value)),
  by = list(vaccination_date = date, geography_code, region)
]

vacc_all <- rbindlist(list(
  vacc_national[, list(
    level = "national", vaccination_date, vaccinated,
    geography = "England", lower_age_limit = NA_integer_
  )],
  vacc_regional[, list(
    level = "regional", vaccination_date, vaccinated,
    geography, lower_age_limit = NA_integer_
  )],
  vacc_local[, list(
    level = "local", vaccination_date, vaccinated,
    geography = geography_code, region, lower_age_limit = NA_integer_
  )],
  vacc_age[, list(
    level = "age_school", vaccination_date, vaccinated,
    geography = "England", lower_age_limit
  )]
), fill = TRUE)
## campaign metrics report a tail of near-zero doses beyond the end of the
## archive, where first, second and third doses are no longer available, so the
## series stops where all sources are present
vacc_all <- vacc_all[vaccination_date <= max(national$date)]

vacc_all[, vaccinated := round(vaccinated)]
setkey(vacc_all, level, geography, lower_age_limit, vaccination_date)
fwrite(vacc_all, here::here("data-processed", "vacc.csv"))
