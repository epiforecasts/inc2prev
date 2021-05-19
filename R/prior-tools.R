#' Tune an inverse Gamma distribtuion to achieve the target truncation
#'
#' @description Allows an inverse gamma distribution to be tuned so that less
#' than 0.01 of its probability mass function falls outside of the specified
#' bounds. This is useful when using an inverse gamma prior, for example for a
#' Gaussian process.
#' @param model A stan model to use to tune the inverse gamma truncation. If
#' not supplied the package default is compiled and used.
#' @param lower Numeric, defaults to 2. Lower truncation bound.
#' @param upper Numeric, defaults to 21. Upper truncation bound.
#'
#' @return A list of alpha and beta values that describe a inverse gamma
#' distribution that achieves the target truncation.
#' @export
#' @importFrom rstan extract sampling stan_model
#' @examples
#' tune_inv_gamma()
tune_inv_gamma <- function(model= inc2prev::load_model("tune_inv_gamma"),
                           lower = 14, upper = 90) {
  # optimise for correct upper and lower probabilites
  fit <- sampling(model,
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

  alpha <- extract(
    fit, "alpha"
  )
  beta <- extract(
    fit, "beta"
  )

  out <- list(
    alpha = round(unlist(unname(alpha)), 1),
    beta = round(unlist(unname(beta)), 1)
  )
  return(out)
}
