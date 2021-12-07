library(dplyr)

ons_to_nhse_region <- function(x) {
  recode(
    x,
    `North East` = "North East and Yorkshire",
    `Yorkshire and The Humber` = "North East and Yorkshire",
    `East Midlands` = "Midlands",
    `West Midlands` = "Midlands",
    )
}

read_cis <- function(fill_missing = TRUE) {
  pops <- read_pop()
  ## get prevalence by ONS region
  prev_regional <- readr::read_csv(here::here("data", "cis.csv")) %>%
    filter(level != "local") %>%
    left_join(pops %>%
              filter(level != "age_school") %>%
              select(geography_code, population),
              by = "geography_code")
  ## get local prevalence, filling with ONS region estimates where missing
  prev_local <- readr::read_csv(here::here("data", "cis.csv")) %>%
    filter(level == "local")
  if (fill_missing) {
    missing_dates_local <- seq(as.Date("2021-03-21"), as.Date("2021-06-27"), by = 7)
    additional_dates <-
      expand_grid(geography_code = unique(prev_local$geography_code),
                  start_date = missing_dates_local) %>%
      mutate(end_date = start_date + 6) %>%
      inner_join(prev_local %>%
                 select(level, geography, geography_code, region) %>%
                 distinct(),
                 by = "geography_code") %>%
      inner_join(prev_regional %>%
                 select(start_date, region = geography,
                        starts_with("proportion")))
    prev_local <- prev_local %>%
      bind_rows(additional_dates) %>%
      arrange(geography_code, start_date)
  }
  ## convert retional prevalence to NHSE regional prevalence
  prev_regional <- prev_regional %>%
  select(level, start_date,
         end_date,
         middle = proportion_pos,
         lower = proportion_pos_low_95,
         upper = proportion_pos_high_95,
         variable = geography,
         population) %>%
    mutate(variable = ons_to_nhse_region(variable)) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(value = sum(population * value) / sum(population),
              population = sum(population),
              .groups = "drop") %>%
    pivot_wider()
  ## finalise local prevalence
  prev_local <- prev_local %>%
    left_join(pops %>%
              filter(level != "age_school") %>%
              select(geography_code, population),
              by = "geography_code") %>%
    select(level,
           start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = geography_code, region,
           population)
  prev_age <- readr::read_csv(here::here("data", "cis_age.csv")) %>%
    left_join(pops %>%
              filter(level == "age_school") %>%
              select(lower_age_limit, population),
              by = "lower_age_limit") %>%
    mutate(age_group = limits_to_agegroups(lower_age_limit)) %>%
    select(level,
           start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = age_group,
           population)
  prev <- bind_rows(prev_regional, prev_local, prev_age) %>%
    mutate(date = start_date + (end_date - start_date) / 2)

  return(prev)
}

read_pop <- function() {
  readr::read_csv(here::here("data", "populations.csv"))
}

read_ab <- function() {
  pops <- read_pop()
  lower_age_limits <- read_cis() %>%
    filter(level == "age_school") %>%
    select(variable) %>%
    distinct() %>%
    mutate(lower_age_limit = parse_number(sub("-\\+.*$", "", variable))) %>%
    pull(lower_age_limit)
  ab_regional <- readr::read_csv(here::here("data", "ab.csv")) %>%
    left_join(pops %>%
              filter(level != "age_school") %>%
              select(geography_code, population),
              by = "geography_code") %>%
    select(level, start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = geography,
           population) %>%
    mutate(variable = ons_to_nhse_region(variable)) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(value = sum(population * value) / sum(population),
              population = sum(population),
              .groups = "drop") %>%
    pivot_wider()
  ab_age <- readr::read_csv(here::here("data", "ab_age.csv")) %>%
    left_join(pops %>%
              filter(level == "age_school") %>%
              select(lower_age_limit, population),
              by = "lower_age_limit") %>%
    mutate(level = "age_school") %>%
    select(level,
           start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = lower_age_limit,
           population) %>%
    mutate(variable =
             reduce_agegroups(variable, lower_age_limits)) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(value = sum(population * value) / sum(population),
              population = sum(population),
              .groups = "drop") %>%
    pivot_wider() %>%
    mutate(variable = limits_to_agegroups(variable))
  ab <- bind_rows(ab_regional, ab_age) %>%
    mutate(date = start_date + (end_date - start_date) / 2)
  return(ab)
}

read_vacc <- function() {
  pops <- read_pop()
  vacc_regional <- readRDS(here::here("data", "vacc.rds")) %>%
    filter(level != "age_school") %>%
    left_join(pops %>%
              filter(level != "age_school") %>%
              select(geography, population),
              by = "geography") %>%
    mutate(vaccinated = vaccinated / population) %>%
    select(level, date = vaccination_date,
           vaccinated, variable = geography)
  vacc_age <- readRDS(here::here("data", "vacc.rds")) %>%
    filter(level == "age_school") %>%
    left_join(pops %>%
              filter(level == "age_school") %>%
              select(lower_age_limit, population),
              by = "lower_age_limit") %>%
    mutate(vaccinated = vaccinated / population,
           age_group = limits_to_agegroups(lower_age_limit)) %>%
    select(level, date = vaccination_date,
           vaccinated, variable = age_group)
  vacc <- bind_rows(vacc_regional, vacc_age)
  return(vacc)
}

read_early <- function() {
  early <- readr::read_csv(here::here("data", "early-seroprevalence.csv"))
  return(early)
}

