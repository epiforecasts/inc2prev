library(data.table)
library(EpiNow2)

# define required stan data
stan_data <- function(prev, prob_detectable, ut = 14, region = "England",
                      gt = list(
                        mean = 3.64, mean_sd = 0.71, sd = 3.08,
                        sd_sd = 0.77, max = 15
                      ),
                      gp_m = 0.3, gp_ls = c(7, 60)) {
  # nolint start
  # extract a single region for prevalence and build features
  prev <- copy(prev)[geography %in% region][, .(date, prev = middle)]
  prev[, time := date - min(date)]

  # summarise prob_detectable for simplicity
  prob_detectable <- melt(
    copy(prob_detectable),
    value.name = "p", id.vars = "sample"
  )
  prob_detectable <- prob_detectable[, .(p = median(p)), by = variable]
  prob_detectable[, time := as.numeric(as.character(variable))]

  # build stan data
  dat <- list(
    ut = ut,
    ot = max(prev$time),
    t = ut + max(prev$time),
    obs = length(prev$prev),
    prev = prev$prev,
    prev_time = prev$time,
    prob_detect = prob_detectable$p,
    pbt = max(prob_detectable$time)
  )

  # gaussian process parameters

  dat$M <- ceiling(dat$t * gp_m)
  dat$L <- 2
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$t
  }
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
