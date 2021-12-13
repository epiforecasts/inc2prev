library(ggplot2)
library(here)
library(data.table)
library(dplyr)
library(truncnorm)
library(forcats)
library(socialmixr)
library(readr)
library(purrr)

## Get tools
functions <- list.files(here("R"), full.names = TRUE)
walk(functions, source)

geo_levels <- c("national", "regional", "local")
other_levels <- c("age_school")

levels <- c(geo_levels, other_levels)

prev <- read_cis()
estimates <- readRDS(here::here("outputs", "estimates.rds"))
samples <- readRDS(here::here("outputs", "samples.rds"))
early <- read_early()

safe_plot_wrapper <- purrr::safely(plot_wrapper)
map(
  levels, plot_wrapper,
  prev = prev, estimates = estimates,
  samples = samples, early = early
)
