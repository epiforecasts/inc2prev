suppressMessages(library(ggplot2))
suppressMessages(library(here))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(truncnorm))
suppressMessages(library(forcats))
suppressMessages(library(socialmixr))
suppressMessages(library(readr))
suppressMessages(library(purrr))
suppressMessages(library(docopt))

doc <- "
Estimate incidence from ONS positivity prevalence data,
possibly including antibody and vaccination data
Usage:
    estimate.R

Options:
    -h --help Show this screen
    -a --ab   Use antibody data
"

## if running interactively can set opts to run with options
if (interactive()) {
  if (!exists("opts")) opts <- list()
} else {
  opts <- docopt(doc)
}

antibodies <- !is.null(opts$ab) && opts$ab

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

geo_levels <- c("national", "regional", "local")
other_levels <- c("age_school")

levels <- c(geo_levels, other_levels)

prev <- read_cis()
if (antibodies) {
  ab <- read_ab()
} else {
  ab <- NULL
}

suffix <- ifelse(antibodies, "_ab", "")
estimates <- readRDS(paste0("outputs/estimates", suffix, ".rds"))
samples <- readRDS(paste0("outputs/samples", suffix, ".rds"))
early <- read_early()

safe_plot_wrapper <- purrr::safely(plot_wrapper)
map(
  levels, plot_wrapper,
  prev = prev, ab = ab, estimates = estimates,
  samples = samples, early = early
)
