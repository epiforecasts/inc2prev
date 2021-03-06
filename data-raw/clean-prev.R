# By Joel Hellewell.
library(data.table)
library(reshape2)
library(readxl)
library(dplyr)
library(lubridate)

# table of interest
table <- "data-raw/20210324_Reference_Table.xlsx"

# Read in data
ons_dt <- readxl::read_xlsx(table,
  sheet = "1i", skip = 6, col_types = c("date", "date", rep("numeric", 27))
)

ons_eng <- readxl::read_xlsx(table,
  sheet = "1d", skip = 7, col_types = c("date", "date", rep("numeric", 13))
)[, 1:5]

colnames(ons_eng) <- c("start_date", "end_date", "middle", "lower", "upper")

ons_eng$geography <- "England"

ons_eng <- as.data.table(ons_eng)
ons_eng <- ons_eng[!is.na(start_date)]

places <- c(
  "North East", "North West", "Yorkshire and the Humber",
  "East Midlands", "West Midlands", "East of England",
  "London", "South East", "South West"
)

colnames(ons_dt) <- c(
  "start_date", "end_date",
  paste(
    rep(c("midle", "lower", "upper"), length(places)),
    rep(places, rep(3, length(places)))
  )
)

ons_dt <- reshape2::melt(ons_dt, id = c("start_date", "end_date"))

ons_dt <- ons_dt %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(
    geography = substr(variable, start = 7, stop = length(variable)),
    type = substr(variable, 1, 5)
  )


ons_dt$variable <- NULL

ons_dt <- reshape2::dcast(
  ons_dt,
  formula = start_date + end_date + geography ~ type
)

colnames(ons_dt)[5] <- "middle"

ons_dt$geography <- as.factor(ons_dt$geography)

ons_dt <- as.data.table(ons_dt)

# Join England to regions
ons_dt <- data.table::rbindlist(list(ons_eng, ons_dt), use.names = TRUE)

# Variable for midpoint of 14 day intervals
ons_dt[, date := lubridate::ymd(ons_dt$start_date) + 7]

# save
data.table::fwrite(ons_dt, "data/ons-prev.csv")
