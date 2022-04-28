extract_publication_dates <- function(x) {
  ## extract publication dates
  publication_dates <- sub("^.+([0-9]{8}?)[^/]+$", "\\1", x)
  ## manually correct misnamed files
  publication_dates[grep("v3", publication_dates)] <- "20200522"
  publication_dates[grep("v2", publication_dates)] <- "20200521"
  publication_dates[grep("v1", publication_dates)] <- "20200514"
  swapped_dates <- which(as.integer(substr(publication_dates, 5, 6)) > 12)
  publication_dates[swapped_dates] <- paste0(
    substr(publication_dates[swapped_dates], 5, 8),
    substr(publication_dates[swapped_dates], 3, 4),
    substr(publication_dates[swapped_dates], 1, 2)
  )
  ## fix typo
  swapped_year <- grepl("^2202", publication_dates)
  publication_dates[swapped_year] <- paste0(
    "2022", substr(publication_dates[swapped_year], 5, 8)
  )
  publication_dates <- ymd(publication_dates)
}
