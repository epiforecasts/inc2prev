library(dplyr)

ons_to_nhse_region <- function(x) {
  recode(
    x,
    `North East` = "North East and Yorkshire",
    `Yorkshire and the Humber` = "North East and Yorkshire",
    `East Midlands` = "Midlands",
    `West Midlands` = "Midlands",
  )
}

read_cis <- function() {
  pops <- read_pop()
  prev_regional <- readr::read_csv(here::here("data", "cis.csv")) %>%
    filter(level != "local") %>%
    left_join(pops %>%
      filter(level != "age_school") %>%
      select(geography_code, population),
    by = "geography_code"
    ) %>%
    select(level, start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = geography,
      population
    ) %>%
    mutate(variable = ons_to_nhse_region(variable)) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(
      value = sum(population * value) / sum(population),
      population = sum(population),
      .groups = "drop"
    ) %>%
    pivot_wider()
  prev_local <- readr::read_csv(here::here("data", "cis.csv")) %>%
    filter(level == "local") %>%
    left_join(pops %>%
      filter(level != "age_school") %>%
      select(geography_code, population),
    by = "geography_code"
    ) %>%
    select(level,
      start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = geography_code, region,
      population
    )
  prev_age <- readr::read_csv(here::here("data", "cis_age.csv")) %>%
    left_join(pops %>%
      filter(level == "age_school") %>%
      select(lower_age_limit, population),
    by = "lower_age_limit"
    ) %>%
    mutate(age_group = limits_to_agegroups(lower_age_limit)) %>%
    select(level,
      start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = age_group,
      population
    )
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
    by = "geography_code"
    ) %>%
    select(level, start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = geography,
      population
    ) %>%
    mutate(variable = ons_to_nhse_region(variable)) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(
      value = sum(population * value) / sum(population),
      population = sum(population),
      .groups = "drop"
    ) %>%
    pivot_wider()
  ab_age <- readr::read_csv(here::here("data", "ab_age.csv")) %>%
    left_join(pops %>%
      filter(level == "age_school") %>%
      select(lower_age_limit, population),
    by = "lower_age_limit"
    ) %>%
    mutate(level = "age_school") %>%
    select(level,
      start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = lower_age_limit,
      population
    ) %>%
    mutate(
      variable =
        reduce_agegroups(variable, lower_age_limits)
    ) %>%
    pivot_longer(c(middle, lower, upper)) %>%
    group_by(level, start_date, end_date, variable, name) %>%
    summarise(
      value = sum(population * value) / sum(population),
      population = sum(population),
      .groups = "drop"
    ) %>%
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
    by = "geography"
    ) %>%
    mutate(vaccinated = vaccinated / population) %>%
    select(level,
      date = vaccination_date,
      vaccinated, variable = geography
    )
  vacc_age <- readRDS(here::here("data", "vacc.rds")) %>%
    filter(level == "age_school") %>%
    left_join(pops %>%
      filter(level == "age_school") %>%
      select(lower_age_limit, population),
    by = "lower_age_limit"
    ) %>%
    mutate(
      vaccinated = vaccinated / population,
      age_group = limits_to_agegroups(lower_age_limit)
    ) %>%
    select(level,
      date = vaccination_date,
      vaccinated, variable = age_group
    )
  vacc <- bind_rows(vacc_regional, vacc_age)
  return(vacc)
}

read_early <- function() {
  early <- readr::read_csv(here::here("data", "early-seroprevalence.csv"))
  return(early)
}

read_prob_detectable <- function() {
  data.table::fread("data/prob_detectable.csv")
}
