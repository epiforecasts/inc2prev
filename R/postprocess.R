i2p_draws <- function(fit, variables = NULL, samples = 100) {
  fit$draws(variables) %>%
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(sample = 1:dplyr::n()) %>%
    dplyr::filter(sample <= samples) %>%
    dplyr::select(-.chain, -.iteration, -.draw) %>%
    tidyr::pivot_longer(matches("[0-9]")) %>%
    dplyr::mutate(
      index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name)
    )
}

i2p_summarise <- function(fit, variables = NULL,
                          quantiles = seq(0.05, 0.95, by = 0.05)) {
  fit$summary(
    variables = variables, ~ quantile(.x, probs = quantiles)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name)
    )
}

## translate index into date
i2p_add_dates <- function(dt, prev, data) {
  start_date <- min(prev$start_date, na.rm = TRUE)
  dates <- prev$date
  ut <- data$ut

  dt <- data.table::as.data.table(dt)[
    ,
    date := fcase(
      name %in% c("infections", "dcases", "dab", "pop_prev"),
      index - 1 + start_date - ut,
      name == "est_prev", dates[index],
      name == "r", index + start_date - ut,
      name == "R", index - 1 + start_date
    )
  ]
  return(dt[])
}
