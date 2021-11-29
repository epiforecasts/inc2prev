read_cis <- function() {
  prev_regional <- readRDS(here::here("data", "cis.rds")) %>%
    filter(level != "local") %>%
    select(level, start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = geography,
           population)
  prev_local <- readRDS(here::here("data", "cis.rds")) %>%
    filter(level == "local") %>%
    select(level,
           start_date,
           end_date,
           middle = proportion_pos,
           lower = proportion_pos_low_95,
           upper = proportion_pos_high_95,
           variable = geography_code, region,
           population)
  prev_age <- readRDS(here::here("data", "cis_age.rds")) %>%
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
