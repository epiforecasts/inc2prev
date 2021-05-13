#' Expose internal package stan functions in R
#'
#' @description This function exposes internal stan functions in R from a user
#' supplied list of target files. Allows for testing of stan functions in R and
#' potentially user use in R code.
#' @source https://github.com/epiforecasts/EpiNow2/
#' @param files A character vector indicating the target files
#' @param target_dir A character string indicating the target directory for the
#' file
#' @param ... Additional arguments passed to `rstan::expose_stan_functions`.
#' @return NULL
#' @export
#' @importFrom rstan expose_stan_functions stanc
#' @importFrom purrr map_chr
expose_stan_fns <- function(files, target_dir, ...) {
  functions <- paste0(
    "\n functions{ \n",
    paste(purrr::map_chr(
      files,
      ~ paste(readLines(file.path(target_dir, .)), collapse = "\n")
    ),
    collapse = "\n"
    ),
    "\n }"
  )
  expose_stan_functions(stanc(model_code = functions), ...)
  return(invisible(NULL))
}

load_model <- function(model = "inc2prev", ...) {
  model <- match.arg(model, choics = c("inc2prev", "tune_inv_gamma"))
  stan_model(
    file = system.file(paste0("stan/", model, ".stan"), package = "inc2prev"),
    ...
  )
}

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
tune_inv_gamma <- function(model, lower = 14, upper = 90) {
  if (missing(model)) {
      model <- load_model("tune_inv_gamma")
  }
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

  alpha <- extract(fit, "alpha")
  beta <- extract(fit, "beta")

  out <- list(
    alpha = round(unlist(unname(alpha)), 1),
    beta = round(unlist(unname(beta)), 1)
  )
  return(out)
}