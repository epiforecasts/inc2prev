library(data.table)
library(EpiNow2)

# define required stan data
stan_data <- function(prev, prob_detectable, ut = 14, region = "England",
                      population = 56286961,
                      gt = list(
                        mean = 3.64, mean_sd = 0.71, sd = 3.08,
                        sd_sd = 0.77, max = 15
                      ),
                      gp_m = 0.3, gp_ls = c(14, 90)) {
  # nolint start
  # extract a single region for prevalence and build features
  prev <- copy(prev)[geography %in% region]
  prev <- prev[, .(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    date = as.Date(date),
    prev = middle,
    sd = (upper - lower) / (2 * 1.96)
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
    t = ut + max(prev$etime),
    obs = length(prev$prev),
    prev = prev$prev,
    prev_sd2 = prev$sd^2,
    prev_time = prev$time,
    prev_stime = prev$stime,
    prev_etime = prev$etime,
    prob_detect_mean = rev(prob_detectable$mean),
    prob_detect_sd = rev(prob_detectable$sd),
    pbt = max(prob_detectable$time) + 1,
    N = population,
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

plot_prev <- function(fit, prev) {
  fit$summary(
    variables = "est_prev",
    ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))
  ) %>%
    mutate(
      date = prev$date + 3
    ) %>%
    ggplot() +
    aes(x = date, y = `50%`, ymin = `5%`, ymax = `95%`) +
    geom_linerange(
      data = prev, aes(y = NULL, ymin = lower, ymax = upper),
      size = 1.1
    ) +
    geom_point(
      data = prev, aes(y = middle, ymin = NULL, ymax = NULL),
      col = "black", size = 1.1
    ) +
    geom_linerange(col = "lightblue", size = 1.1) +
    geom_point(col = "#0eace0", size = 1.3) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    theme_minimal()
}
