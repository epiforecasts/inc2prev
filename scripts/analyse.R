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
Analyse outputs of the inc2prev model
Usage:
    estimate.R [--ab] [--local | --age | --variants]
    estimate.R -h | --help

Options:
    -h --help        Show this screen
    -a --ab          Use antibody data
    -l, --local      Analyse local dynamics
    -g, --age        Analyse age
    -v, --variants   Analyse variants
"

## if running interactively can set opts to run with options
if (interactive()) {
  if (!exists("opts")) opts <- list()
} else {
  opts <- docopt(doc)
}

antibodies <- !is.null(opts$ab) && opts$ab
local <- !is.null(opts$local) && opts$local
age <- !is.null(opts$age) && opts$age
variants <- !is.null(opts$variants) && opts$variants

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

if (local) {
  suffix <- "_local"
} else if (age) {
  suffix <- "_age"
} else if (variants) {
  suffix <- "_variants"
} else {
  suffix <- "_regional"
}
suffix <- paste0(suffix, ifelse(antibodies, "_ab", ""))

estimates <- readRDS(paste0("outputs/estimates", suffix, ".rds"))
samples <- readRDS(paste0("outputs/samples", suffix, ".rds"))

nhse <- "Midlands" %in% estimates$variable
prev <- read_cis(nhse_regions = nhse)
if (antibodies) {
  ab <- read_ab(nhse_regions = nhse)
} else {
  ab <- NULL
}

if (variants) {
  early <- NULL
} else {
  early <- read_early(nhse_regions = nhse)
}

levels <- unique(estimates$level)

updated_samples <- map(
  levels, plot_wrapper,
  prev = prev, ab = ab, estimates = estimates,
  samples = samples, early = early,
  suffix = suffix, extension = ".svg"
)

cumulative <- bind_rows(updated_samples) %>%
  filter(name == "cumulative_infections")

saveRDS(cumulative,
        here::here("outputs", paste0("cumulative", suffix, ".rds")))

