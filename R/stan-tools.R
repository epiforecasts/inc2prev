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
  files <- purrr::map_chr(files, ~ file.path(target_dir, .))
  functions <- stan_block("functions", files)
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

stan_block <- function(block, contents) {
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
                          generated_quantities = c()) {

  model <- paste(
    stan_block("functions", functions),
    stan_block("data", data),
    stan_block("transformed data", tdata),
    stan_block("parameters", parameters),
    stan_block("transformed parameters", tparameters),
    stan_block("model", model),
    stan_block("generated quantities", generated_quantities),
    collapse = "\n"
  )
  return(model)
}