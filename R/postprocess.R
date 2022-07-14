i2p_draws <- function(fit, variables = NULL, samples = 100) {
  fit$draws(variables) %>%
    posterior::as_draws_df() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(sample = 1:dplyr::n()) %>%
    dplyr::filter(sample <= samples) %>%
    dplyr::select(-.chain, -.iteration, -.draw) %>%
    tidyr::pivot_longer(matches("[0-9]")) %>%
    dplyr::mutate(
      n_index = as.integer(sub("^.*\\[([0-9]+),.+]$", "\\1", name)),
      t_index = as.integer(sub("^.*\\[.+,([0-9]+)]$", "\\1", name)),
      p_index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name)
    )
}

i2p_summarise <- function(fit, variables = NULL,
                          quantiles = seq(0.05, 0.95, by = 0.05)) {
  fit$summary(
    variables = variables, mean, sd, median, mad, ~ posterior::quantile2(.x, probs = quantiles)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(name = variable) %>%
    dplyr::mutate(
      n_index = as.integer(sub("^.*\\[([0-9]+),.+]$", "\\1", name)),
      t_index = as.integer(sub("^.*\\[.+,([0-9]+)]$", "\\1", name)),
      p_index = as.integer(sub("^.*\\[([0-9]+)]$", "\\1", name)),
      name = sub("\\[.*$", "", name)
    )
}

## translate index into date
i2p_add_date <- function(dt, prev, ab, data) {
  start_date <- min(prev$start_date, na.rm = TRUE)
  ut <- data$ut
  if (is.null(ab)) ab <- prev ## needed to avoid error from fcase below

  dt <- suppressWarnings(
    data.table::as.data.table(dt)[
      ,
      date := fcase(
        name %in% c("infections", "dcases", "dab", "gen_dab"),
        t_index - 1 + start_date - ut,
        name == "est_prev", prev$date[t_index],
        name == "est_ab", ab$date[t_index],
        name == "r", t_index - 1 + start_date,
        name == "R", t_index - 1 + start_date
      )
    ]
  )
  return(dt[])
}

## translate index into variable
i2p_add_var <- function(dt, prev, data, var_col = NULL) {
  if (is.null(var_col)) return(dt[])
  vars <- unique(rownames(data$prev))
  ab_index <- data$ab_index

  dt <- suppressWarnings(
    data.table::as.data.table(dt)[
      ,
      paste(var_col) := fcase(
        name %in% c("est_prev", "est_ab", "infections",
		    "dcases", "dab", "gen_dab", "r", "R",
		    "eta", "init_growth"),
	vars[n_index],
        name %in% c("rho", "alpha", "init_inc", "init_dab"),
        vars[p_index],
        name %in% c("beta", "gamma", "delta", "k", "l"),
        paste0(vars[ab_index], collapse = ";")
      )
    ]
  )
  return(dt[])
}
