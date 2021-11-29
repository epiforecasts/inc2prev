#' Tune an Inverse Gamma to Achieve the Target Truncation
#'
#' @description `r lifecycle::badge("questioning")`
#' Allows an inverse gamma distribution to be. tuned so that less than 0.01 of its
#' probability mass function falls outside of the specified
#' bounds. This is required when using an inverse gamma prior, for example for a
#' Gaussian process. As no inverse gamma priors are currently in use and this function
#' has some stability issues it may be deprecated at a later date.
#' @param lower Numeric, defaults to 2. Lower truncation bound.
#' @param upper Numeric, defaults to 21. Upper truncation bound.
#'
#' @return A list of alpha and beta values that describe a inverse gamma
#' distribution that achieves the target truncation.
#' @export
#'
#' @examples
#'
#' tune_inv_gamma(lower = 2, upper = 21)
tune_inv_gamma <- function(lower = 2, upper = 21) {
  model <- rstan::stan_model("stan/tune_inv_gamma.stan")
  # optimise for correct upper and lower probabilites
  fit <- rstan::sampling(model,
    data = list(
      u = upper,
      l = lower
    ),
    iter = 1,
    warmup = 0,
    chains = 1,
    algorithm = "Fixed_param",
    refresh = 0
  )

  alpha <- rstan::extract(fit, "alpha")
  beta <- rstan::extract(fit, "beta")

  out <- list(
    alpha = round(unlist(unname(alpha)), 1),
    beta = round(unlist(unname(beta)), 1)
  )
  return(out)
}