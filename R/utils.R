library(data.table)
library(EpiNow2)

# define required stan data
stan_data <- function(prev, prob_detectable, ut = 14,
                      population = 56286961,
                      gt = list(
                        mean = 3.64, mean_sd = 0.71, sd = 3.08,
                        sd_sd = 0.77, max = 15
                      ),
                      gp_m = 0.3, gp_ls = c(14, 90)) {
  # extract a prevalence and build features
  prev <- data.table(prev)[, .(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    date = as.Date(date),
    prev = middle,
    sd = (upper - lower) / (2 * 1.96),
    population
  )]
  prev[, `:=`(
    time = as.integer(date - min(start_date)),
    stime = as.integer(start_date - min(start_date)),
    etime = as.integer(end_date - min(start_date))
  )]

  # summarise prob_detectable for simplicity
  prob_detectable <- melt(
    copy(prob_detectable),
    value.name = "p", id.vars = "sample"
  )
  prob_detectable[, time := as.numeric(as.character(variable))]
  prob_detectable <- prob_detectable[, .(
    median = median(p),
    mean = mean(p),
    sd = sd(p)
  ),
  by = time
  ]
  prob_detectable <- prob_detectable[,
    purrr::map(.SD, signif, digits = 3),
    .SDcols = c("mean", "median", "sd"),
    by = time
  ]
  # define baseline incidence
  baseline_inc <- prev$prev[1] * prob_detectable$mean[ut]

  # build stan data
  dat <- list(
    ut = ut,
    ot = max(prev$etime),
    t = ut + length(prev$prev),
    obs = length(prev$prev),
    prev = prev$prev,
    prev_sd2 = prev$sd^2,
    prev_time = prev$time,
    prev_stime = prev$stime,
    prev_etime = prev$etime,
    prob_detect_mean = rev(prob_detectable$mean),
    prob_detect_sd = rev(prob_detectable$sd),
    pbt = max(prob_detectable$time) + 1,
    N = unique(prev$population),
    inc_zero = log(baseline_inc / (baseline_inc + 1))
  )

  # gaussian process parameters
  dat$M <- ceiling(dat$t * gp_m)
  dat$L <- 2
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$t
  }
  gp_ls
  lsp <- EpiNow2::tune_inv_gamma(gp_ls[1], gp_ls[2])
  dat$lengthscale_alpha <- lsp$alpha
  dat$lengthscale_beta <- lsp$beta

  # define generation time
  dat$gtm <- unlist(gt[c("mean", "mean_sd")])
  dat$gtsd <- unlist(gt[c("sd", "sd_sd")])
  dat$gtmax <- unlist(gt[c("max")])
  # nolint end
  return(dat)
}

library(truncnorm)
library(purrr)

stan_inits <- function(dat) {
  inits <- function() {
    list(
      eta = array(rnorm(dat$M, mean = 0, sd = 0.1)),
      alpha = array(truncnorm::rtruncnorm(1, mean = 0, sd = 0.1, a = 0)),
      sigma = array(truncnorm::rtruncnorm(1, mean = 0.005, sd = 0.0025, a = 0)),
      rho = array(truncnorm::rtruncnorm(1, mean = 36, sd = 21, a = 14, b = 90)),
      prob_detect = purrr::map2_dbl(
        dat$prob_detect_mean, dat$prob_detect_sd / 10,
        ~ truncnorm::rtruncnorm(1, a = 0, b = 1, mean = .x, sd = .y)
      )
    )
  }
  return(inits)
}

library(dplyr)
library(ggplot2)

# plot trend with date over time
plot_trend <- function(fit, var, date_start) {
  fit$summary(
    variables = var,
    ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))
  ) %>%
    mutate(
      time = 1:n(),
      date = date_start + time - 1
    ) %>%
    ggplot() +
    aes(x = date, y = `50%`, ymin = `5%`, ymax = `95%`) +
    geom_line(col = "lightblue", size = 1.4) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = "lightblue", size = 0.6
    ) +
    geom_ribbon(
      fill = "lightblue", alpha = 0.4,
      col = NA, aes(ymin = `20%`, ymax = `80%`)
    ) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    theme_minimal()
}

plot_trace <- function(samples, var, alpha = 0.05) {

  long_samples <- samples %>%
    filter(name == var) %>%
    pivot_longer(matches("^[0-9]+$"), names_to = "sample")

  plot <- long_samples %>%
    ggplot() +
    aes(x = date, y = value, group = sample) +
    geom_line(alpha = alpha) +
    theme_minimal() +
    labs(x = "Date") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %d") +
    facet_wrap(~ variable)

  return(plot)
}


library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

plot_prev <- function(estimates, samples, data, alpha = 0.05,
                      data_source = "ONS Prevalence") {
  trace_plot <- plot_trace(
    samples,
    "pop_prev",
    alpha = alpha
  )

  summary_prev <- estimates %>%
    filter(name == "est_prev") %>%
    mutate(
      middle = `50%`,
      lower = `5%`,
      upper = `95%`,
      type = "Modelled"
    ) %>%
    bind_rows(data %>%
      mutate(
        type = "Estimate"
      ))

  trace_plot +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Prevalence", x = "Date") +
    geom_linerange(
      data = summary_prev,
      aes(
        y = NULL, ymin = lower, ymax = upper, group = NULL,
        col = type
      ),
      size = 1, alpha = 0.2
    ) +
    geom_point(
      data = summary_prev,
      aes(
        y = middle, ymin = NULL, ymax = NULL, group = NULL,
        col = type
      ), size = 1.1, alpha = 0.2
    ) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Dark2") +
    guides(col = guide_legend(title = data_source))
}

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

plot_ltla <- function(estimates, areas, names = c(), var = "pop_prev", days = 60, var_name = "Prevalence") {
  estimates <- estimates %>%
    filter(name == {{var}}) %>%
    filter(date > max(date) - {{days}})
  if (length(names) > 0) {
    search_str <- paste0("(",
                         paste(names, collapse = "|"),
                         ")")
    areas <- areas %>%
      mutate(highlighted = grepl(search_str, ltla_name)) %>%
      group_by(geography_code) %>%
      summarise(highlighted = any(highlighted), .groups = "drop") %>%
      mutate(highlighted = if_else(highlighted, "yes", "no"),
             highlighted = factor(highlighted, levels = c("yes", "no")))
  }
  estimates <- estimates %>%
    inner_join(areas %>% rename(variable = geography_code), by = "variable")
  aesthetics <- list(x = "date",
                     y = "`50%`",
                     group = "variable")
  if (length(names) > 0) {
    aesthetics[["colour"]] <- "highlighted"
    aesthetics[["alpha"]] <- "highlighted"
  } else {
    aesthetics[["colour"]] <- "region"
  }
  p <- ggplot(estimates, mapping = do.call(aes_string, aesthetics)) +
    geom_line()
  if (length(names) > 0) {
    p <- p +
      scale_colour_manual("", values = c("red", "black")) +
      scale_alpha_manual("", values = c(1, 0.25)) +
      theme(legend.position = "none")
  } else {
    p <- p +
      scale_colour_brewer("Region", palette = "Paired")
  }
  p <- p +
    theme_minimal() +
    xlab("") +
    ylab(var_name)
  return(p)
}
