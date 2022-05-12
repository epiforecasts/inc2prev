#' @export
i2p_simulate <- function(dat, param, nsamples = NULL, type = "estimate") {

  type <- match.arg(type, choices = c("estimate", "observations"))
  files <- c("ab.stan", "gaussian_process.stan",
             "observed_in_window.stan", "convolve.stan")

  functions <-
    paste0("\n functions{ \n", paste(purrr::map_chr(
      files, ~paste(readLines(here::here("stan", "functions", .)),
                    collapse = "\n")),
      collapse = "\n"), "\n }")

  ## fix rstan / cmdstanr compatibility issue
  functions <- gsub("array\\[\\] ([a-z]+) ([a-z_]+)", "\\1[] \\2", functions)

  rstan::expose_stan_functions(rstan::stanc(model_code = functions))

  ## transformed data

  PHI <- setup_gp(dat$M, dat$L, dat$t - dat$diff_order)
  vacc_with_ab <- t(apply(dat$vacc, 1, convolve, dat$vacc_ab_delay))

  gp <- lapply(seq_len(dat$n), function(i) {
    upd <- update_gp(PHI, dat$M, dat$L, param$alpha[i], param$rho[i], 
		     param$eta[i, ], 0)
    return(upd)
  })
  gp <- do.call(rbind, gp)

  infections <- lapply(seq_len(dat$n), function(i) {
    inf <- plogis(param$init_inc[i] + gp[i, ])
    return(inf)
  })
  infections <- do.call(rbind, infections)

  dcases <- lapply(seq_len(dat$n), function(i) {
    dca <- convolve(infections[i, ], param$prob_detect)
    return(dca)
  })
  dcases <- do.call(rbind, dcases)

  odcases <- lapply(seq_len(dat$n), function(i) {
    odc <- observed_in_window(dcases[i, ], dat$prev_stime, dat$prev_etime, 
			      dat$ut, dat$obs)
    return(odc)
  })
  odcases <- do.call(rbind, odcases)

  infs_with_potential_abs <- lapply(seq_len(dat$n_ab), function(i) {
    inf <- convolve(infections[i, ], dat$inf_ab_delay)
    return(inf)
  })
  infs_with_potential_abs <- do.call(rbind, infs_with_potential_abs)

  dab <- lapply(seq_len(dat$n), function(i) {
    da <- detectable_antibodies(infs_with_potential_abs[i, ],
  			        vacc_with_ab[i, ],
			        param$beta, param$gamma, param$delta,
			        param$k, param$l, 
			        param$init_dab[param$ab_index[i]], dat$t)
    return(da)
  })
  dab <- do.call(rbind, dab)

  if (type == "observation") {
    odab <- lapply(seq_len(dat$n), function(i) {
      oda <- observed_in_window(dab[i, ], dat$ab_stime, dat$ab_etime, dat$ut, dat$ab_obs)
      return(oda)
    })
    odab <- do.call(rbind, odab)

    combined_ab_sigma <- lapply(seq_len(dat$n), function(i) {
      cs <- sqrt(rep(param$ab_sigma^2, ncol(dat$ab_sd2)) + dat$ab_sd2[i, ])
      return(cs)
    })
    combined_ab_sigma <- do.call(rbind, combined_ab_sigma)

    ab <- rnorm(length(odab), odab, combined_ab_sigma)

    return(ab)
  } else {
    return(dab)
  }
}
