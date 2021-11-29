read_cis <- function() {
  prev_regional <- readr::read_csv(here::here("data", "cis.rds")) %>%
    dplyr::filter(level != "local") %>%
    dplyr::select(level, start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = geography,
      population
    )
  prev_local <- readr::read_csv(here::here("data", "cis.rds")) %>%
    dplyr::filter(level == "local") %>%
    dplyr::select(level,
      start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = geography_code, region,
      population
    )
  prev_age <- readr::read_csv(here::here("data", "cis_age.rds")) %>%
    dplyr::mutate(
      age_group = socialmixr::limits_to_agegroups(lower_age_limit)
    ) %>%
    dplyr::select(level,
      start_date,
      end_date,
      middle = proportion_pos,
      lower = proportion_pos_low_95,
      upper = proportion_pos_high_95,
      variable = age_group,
      population
    )
  prev <- dplyr::bind_rows(prev_regional, prev_local, prev_age) %>%
    dplyr::mutate(date = start_date + (end_date - start_date) / 2)

  return(prev)
}
