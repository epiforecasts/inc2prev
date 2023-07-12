#' @title Expose Stan functions
#' @description Exposes stan functions to R
#' @param files Stan files containing functions to expose
#' @param target_dir Directory containing Stan files
#' @param ... Additional arguments passed to
#'   \code{\link{rstan::expose_stan_functions}}
#' @return NULL (invisible)
#' @family utils
#' @export
#' @importFrom purrr map_chr
#' @importFrom rstan expose_stan_functions stanc
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
  rstan::expose_stan_functions(rstan::stanc(model_code = functions), ...)
  return(invisible(NULL))
}

#' @title Inverse logit function
#' @description The inverse logit function
#' @param x a logit values or vector of logit values
#' @return the probability calculated as inverse logit
#' @family utils
#' @examples
#' inv_logit(c(-10, 1, 0, 100))
#' @export
inv_logit <- function(x) {
  il <- 1 / (1 + exp(-x))
  return(il)
}

#' @title Logit function
#' @description The logit function
#' @param p a probability or vector of probabilities
#' @return the logit of the probability
#' @family utils
#' @examples
#' logit(c(0.01, 0.1, 0.5, 0.9, 1))
#' @export
logit <- function(p) {
  l <- log(p / (1 - p))
  return(l)
}
