library(boot)
library(data.table)
library(purrr)

dt <- fread("data-raw/hellewell-prI-posterior.csv")

prob_cont <- function(time_since_inf) {
  . <- NULL
  p <- NULL
  prob_detect <- copy(dt)[
    ,
    p := boot::inv.logit(
      beta1 + beta2 * (time_since_inf - cutpoint) +
        (time_since_inf - cutpoint) * beta3 * beta2 * fifelse(
          time_since_inf - cutpoint > 0, 1, 0
        )
    )
  ]
  prob_detect <- prob_detect[
    ,
    .(sample = 1:.N, time = time_since_inf, p)
  ]
  return(prob_detect)
}

# Can pass a vector of quantiles
prob_detectable <- map_df(0:60, prob_cont)
prob_detectable <- dcast(prob_detectable,
  sample ~ time,
  value.var = "p"
)

usethis::use_data(prob_detectable, overwrite = TRUE)
