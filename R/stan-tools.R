stan_block <- function(block, contents, path = NULL) {
  if (!is.null(path)) {
    contents <- purrr::map_chr(contents, ~ file.path(path, .))
  }
  if (!missing(contents)) {
    if (length(contents) > 0) {
        block <- paste0(
          "\n ", block, "{ \n",
          paste(purrr::map_chr(
            contents,
            ~ paste(readLines(.), collapse = "\n")
            ),
            collapse = "\n"
            ),
            "\n }"
            )
    }
  }else{
    block <- paste0("\n ", block, "{}")
  }

  return(block)
}

stan_template <- function(functions = c(), data = c(), tdata = c(),
                          parameters = c(), tparameters = c(), model = c(),
                          generated_quantities = c(), path = NULL) {

  model <- paste(
    stan_block("functions", functions, path),
    stan_block("data", data, path),
    stan_block("transformed data", tdata, path),
    stan_block("parameters", parameters, path),
    stan_block("transformed parameters", tparameters, path),
    stan_block("model", model, path),
    stan_block("generated quantities", generated_quantities, path),
    collapse = "\n"
  )
  return(model)
}

inc2prev <- function(path = system.file("stan", package = "inc2prev")) {

  model <- stan_template(
    functions = c(
      "functions/gaussian_process.stan",
      "functions/rt.stan",
      "functions/prev.stan",
      "functions/generated_quantities.stan"
    ),
    data = c(
      "data/observations.stan",
      "data/prev_obs_model.stan",
      "data/observation_generation.stan",
      "data/summary_measures.stan",
      "data/gaussian_process.stan"
    ),
    tdata = c(
      "tdata/gaussian_process.stan"
    ),
    parameters = c(
      "parameters/gaussian_process.stan",
      "parameters/prev_obs_model.stan"
    ),
    tparameters = c(
      "tparameters-var-def/gaussian_process.stan",
      "tparameters-var-def/observation_generation.stan",
      "tparameters-var-def/prev_obs_model.stan",
      "tparameters/gaussian_process.stan",
      "tparameters/observation_generation.stan",
      "tparameters/prev_obs_model.stan"
    ),
    model = c(
      "model/gaussian_process.stan",
      "model/prev_obs_model.stan"
    ),
    generated_quantities = c(
      "generated-quantities-var-def/summary_measures.stan",
      "generated-quantities/summary_measures.stan"
    ),
    path = path
  )
}
#' Expose internal package stan functions in R
#'
#' @description This function exposes internal stan functions in R from a user
#' supplied list of target files. Allows for testing of stan functions in R and
#' potentially user use in R code.
#' @source https://github.com/epiforecasts/EpiNow2/
#' @param files A character vector indicating the target files
#' @param path A character string indicating the target directory for the
#' file
#' @param ... Additional arguments passed to `rstan::expose_stan_functions`.
#' @return NULL
#' @export
#' @importFrom rstan expose_stan_functions stanc
#' @importFrom purrr map_chr
expose_stan_fns <- function(files,
                            path = system.file("stan", package = "inc2prev"),
                            ...) {
  functions <- stan_block("functions", files, path)
  expose_stan_functions(
    stanc(
      model_code = functions
    ),
    ...
  )
  return(invisible(NULL))
}

load_model <- function(model = "tune_inv_gamma", ...) {
  model <- match.arg(model, choices = c("tune_inv_gamma"))
  stan_model(
    file = system.file(paste0("stan/", model, ".stan"), package = "inc2prev"),
    ...
  )
}
