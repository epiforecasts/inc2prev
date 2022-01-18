prob_detectable <- function(time_since_inf, dt) {
  . <- NULL
  p <- NULL
  prob_detect <- data.table::copy(dt)[
    ,
    p := boot::inv.logit(
      beta1 + beta2 * (time_since_inf - cutpoint) +
        (time_since_inf - cutpoint + 0.5) * beta3 * beta2 * fifelse(
          time_since_inf - cutpoint + 0.5 > 0, 1, 0
        )
    )
  ]
  prob_detect <- prob_detect[
    ,
    .(sample = 1:.N, time = time_since_inf, p)
  ]
  return(prob_detect)
}