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
suppressMessages(library(inc2prev))

source(here::here("scripts", "read.R"))

doc <- "
Analyse outputs of the inc2prev model
Usage:
    estimate.R [--ab] [--higher] [--local | --regional | --age | --variants]
    estimate.R -h | --help

Options:
    -h --help        Show this screen
    -a --ab          Use antibody data
    -i, --higher     Use higher antibody threshold
    -r, --regional   Analyse regional dynamics
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
higher <- !is.null(opts$higher) && opts$higher
regional <- !is.null(opts$regional) && opts$regional
local <- !is.null(opts$local) && opts$local
age <- !is.null(opts$age) && opts$age
variants <- !is.null(opts$variants) && opts$variants

if (local) {
  suffix <- "_local"
} else if (regional) {
  suffix <- "_regional"
} else if (age) {
  suffix <- "_age"
} else if (variants) {
  suffix <- "_variants"
} else {
  suffix <- "_national"
}

if (antibodies) {
  threshold <- ifelse(higher, "higher", "standard")
  ab <- read_ab(nhse_regions = nhse, threshold = threshold)
   suffix <- paste0(suffix, "_ab")
  if (higher) {
    suffix <- paste0(suffix, "_higher")
  }
} else {
  ab <- NULL
}

estimates <- readRDS(paste0("outputs/estimates", suffix, ".rds"))
samples <- readRDS(paste0("outputs/samples", suffix, ".rds"))

nhse <- "Midlands" %in% estimates$variable
prev <- read_cis(nhse_regions = nhse)

if (variants || local) {
  early <- NULL
} else {
  early <- read_early(nhse_regions = nhse)
}

levels <- unique(estimates$level)

dir.create(here::here("figures/additional"), 
	   showWarnings = FALSE, recursive = TRUE)

updated_samples <- map(
  levels, plot_wrapper,
  prev = prev, ab = ab, estimates = estimates,
  samples = samples, early = early,
  suffix = suffix, extension = ".svg"
)

cumulative <- bind_rows(updated_samples) %>%
  filter(name %in% c("cumulative_infections",
                     "cumulative_exposure"))
saveRDS(cumulative,
        here::here("outputs", paste0("cumulative", suffix, ".rds")))

csum <- cumulative %>%
  group_by(name, date, variable, level) %>%
  summarise(x = quantile(value, seq(0.05, 0.95, by = 0.05)),
	    q = paste0("q", seq(5, 95, by = 5)),
	    .groups = "drop") %>%
  pivot_wider(names_from = "q", values_from = "x")

saveRDS(csum,
        here::here("outputs", paste0("cumulative", suffix, ".csv")))


