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
    estimate.R [--ab] [--higher] [--local | --regional | --age | --variants] [--max-report-date=<date>] 
    estimate.R -h | --help

Options:
    -h --help        Show this screen
    -a --ab          Use antibody data
    -i, --higher     Use higher antibody threshold
    -r, --regional   Analyse regional dynamics
    -l, --local      Analyse local dynamics
    -g, --age        Analyse age
    -v, --variants   Analyse variants
    -m, --max-report-date=<date> Latest report date to use for estimation
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
report_date <- opts$max_report_date
if (!is.null(report_date)) report_date <- as.Date(report_date)

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
  suffix <- paste0(suffix, "_ab")
  threshold <- ifelse(higher, "higher", "standard")
  if (higher) {
    suffix <- paste0(suffix, "_higher")
  }
}

if (!is.null(report_date)) {
  suffix <- paste0(suffix, "_", report_date)
}

estimates <- readRDS(paste0("outputs/estimates", suffix, ".rds"))
samples <- readRDS(paste0("outputs/samples", suffix, ".rds"))

nhse <- "Midlands" %in% estimates$variable
prev <- read_cis(nhse_regions = nhse,
                 max_publication_date = report_date)

if (antibodies) {
  ab <- read_ab(nhse_regions = nhse, threshold = threshold,
		max_publication_date = report_date)
} else {
  ab <- NULL
}

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

write_csv(csum,
          here::here("outputs", paste0("cumulative", suffix, ".csv")))


