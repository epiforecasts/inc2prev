library("here")
library("readxl")
library("rvest")
library("tidyr")
library("dplyr")
library("lubridate")
library("ggplot2")
library("purrr")
library("janitor")
library("socialmixr")
library("readr")

## create directory for antibody data if it doesn't exist
ab_dir <- here::here("data", "ab")
dir.create(ab_dir, showWarnings = FALSE, recursive = TRUE)

## creata URLs that list spreadsheets
years <- c(2021)
urls <- paste0("https://www.ons.gov.uk/peoplepopulationandcommunity/",
               "healthandsocialcare/conditionsanddiseases/datasets/",
               "/coronaviruscovid19antibodydatafortheuk/", years)

## get URLs of the spreadsheets, scraped from the web pages
file_urls <- lapply(urls, function(url) {
  session <- session(url)
  file_url <- session %>%
    html_nodes(xpath = paste0(
                 "//a[contains(concat(' ', ",
                 "normalize-space(@class),' '),' btn--primary ')]"
               )) %>%
    html_attr("href")
  return(file_url)
}) %>%
  unlist() %>%
  grep("\\.xlsx?$", value = TRUE, .)

## construct tibble with files to download
df_dl <- tibble(file_url = file_urls) %>%
  mutate(file_name = sub("^.*/([^/]+)$", "\\1", file_url),
         file_path = file.path(ab_dir, file_name),
         full_url = paste0("https://www.ons.gov.uk", file_url)) %>%
  filter(!file.exists(file_path))

## define levels to extract and table structure
columns <- c(national = 5, regional = 6, age_school = 6)
super_headers <- c(regional = "region", age_school = "lower_age_limit")

## list all files
files <- list.files(here::here("data", "ab"), full.names = TRUE)
list_file <- here::here("data", "ab_files.rds")

## if no new URLs there is nothing to do
if (nrow(df_dl) > 0) {
  df_dl %>%
    rowwise() %>%
    mutate(ret = download.file(full_url, file_path))
  if (any(df_dl$ret != 0)) warning("Some downloads failed")
}

## define geography codes not in data
geography_codes <- c(England = "E92000001",
                     `North East` = "E12000001",
                     `North West` = "E12000002",
                     `Yorkshire and The Humber` = "E12000003",
                     `East Midlands` = "E12000004",
                     `West Midlands` = "E12000005",
                     `East of England` = "E12000006",
                     `London` = "E12000007",
                     `South East` = "E12000008",
                     `South West` = "E12000009")

if (file.exists(list_file) && setequal(files, readRDS(list_file))) {
##  stop("Nothing new to extract")
}

## construct list of data frames with positivity
ab <- list()
for (level in names(columns)) {
  ab[[level]] <- lapply(files, function(x) {
    ## first,  get table of contents sheet to work out which sheet we want
    contents_sheet <- read_excel(x, sheet = "Contents") %>%
      clean_names()
    if (level == "national") {
      contents_sheet <- contents_sheet %>%
        filter(grepl("Antibodies$", contents)) %>%
        head(n = 1)
    } else if (level == "regional") {
      contents_sheet <- contents_sheet %>%
        filter(grepl("by region", contents)) %>%
        head(n = 1)
    } else if (level == "age_school") {
      contents_sheet <- contents_sheet %>%
        filter(grepl("by age group$", contents)) %>%
        head(n = 1)
    } else {
      stop("Unknown level: ", level)
    }
    ## extract table number
    sheet <- sub("^Table ([^ ]+) ?- .*$", "\\1", contents_sheet$contents)
    if (length(sheet) == 1) {
      ## we found the sheet, now we get a preview so we can work out where in
      ## the sheet the actual table is
      preview <- read_excel(x, sheet = sheet) %>%
        remove_empty("cols") %>%
        clean_names()
      ## get row that contains the headers
      headers_row <- min(which(!is.na(preview$x2)))
      skip <- headers_row + if_else(level %in% c("regional", "age_school"), 1, 0)
      if (level %in% c("regional", "age_school")) {
        headers <- preview[headers_row, 2:ncol(preview)] %>%
          t() %>%
          as_tibble(.name_repair = "minimal") %>%
          rename(header = 1) %>%
          fill(header)
        if (level == "age_school") {
          headers <- headers %>%
            mutate(header = sub("^Age ([0-9]+) *- *([0-9]+).*$", "\\1|\\2", header),
                   header = sub("^Age ([0-9]+)\\+", "\\1|", header)) %>%
            separate(header, c("from", "to"), sep = "\\|")
        }
      }
      ## having figured out where the table is and extracted the date, read the table
      if (is.infinite(skip)) return(NULL) ## couldn't find data
      ## find max row to read
      n_max <- min(which(!grepl("^[0-9]", unlist(preview[(skip + 1):nrow(preview), 1])))) - 1
      data <- read_excel(x, sheet = sheet, skip = skip, n_max = n_max, .name_repair = "minimal") %>%
        remove_empty("cols") %>%
        clean_names()
      colnames(data) <- c("weekly_period", sub("_[0-9]+$", "", colnames(data)[2:ncol(data)]))
      if (level == c("regional")) {
        colnames(data)[2:ncol(data)] <-
          paste(colnames(data)[2:ncol(data)], headers$header, sep = "|")
      } else if (level == "age_school") {
        colnames(data)[2:ncol(data)] <-
          paste(colnames(data)[2:ncol(data)], headers$from, sep = "|")
      }
      data <- data[, !duplicated(colnames(data))]
      data <- data %>%
        pivot_longer(2:ncol(.))
      if (level %in% names(super_headers)) {
        data <- data %>%
          separate(name, c("name", super_headers[level]), sep = "\\|")
      }
     data <- data %>%
        pivot_wider() %>%
        separate(weekly_period, c("start_date", "end_date"), sep = " to ") %>%
        mutate(start_date = dmy(start_date),
               end_date = dmy(end_date))
      data <- data %>%
        select(1:columns[level]) %>%
        rename(proportion_pos = columns[[level]] - 2,
               proportion_pos_low_95 = columns[[level]] - 1,
               proportion_pos_high_95 = columns[[level]]) %>%
        mutate_at(vars(starts_with("proportion")), as.numeric) %>%
        pivot_longer(starts_with("proportion")) %>%
        replace_na(list(value = 0)) %>%
        mutate(value = value / 100) %>%
        pivot_wider()
      if (level %in% c("national", "regional", "age_school")) {
        if (level %in% c("national", "age_school")) {
          data <- data %>%
            mutate(region = NA_character_,
                   geography = "England")
        } else if (level == "regional") {
          data <- data %>%
            mutate(region = sub("the Humber", "The Humber", region)) %>%
            mutate(geography = region)
        }
        data <- data %>%
          mutate(geography_code = geography_codes[geography])
      }
      if (level == "age_school") {
        data <- data %>%
          mutate(lower_age_limit = as.integer(lower_age_limit)) %>%
          select(start_date, end_date, geography, geography_code,
                 lower_age_limit, starts_with("proportion"))
      } else {
        data <- data %>%
          select(start_date, end_date, geography, geography_code,
                 starts_with("proportion"))
      }
      return(data %>%
             mutate(file_name = x))
    } else {
      return(NULL)
    }
  })
  ab[[level]] <- ab[[level]] %>%
    bind_rows() %>%
    distinct() %>% ## avoid duplicate rows
    mutate(level = level)
}

## combine it all into one data frame
combined <- ab %>%
  bind_rows() %>%
  group_by(file_name) %>%
  mutate(report_date = max(end_date)) %>%
  ungroup() %>%
  filter(report_date == max(report_date)) %>%
  arrange(level, geography, start_date)

## save
write_csv(combined %>%
          filter(level %in% c("national", "regional")) %>%
          remove_empty(which = "cols"),
          here::here("data", "ab.csv"))
write_csv(combined %>%
          filter(level %in% c("age_school")) %>%
          remove_empty(which = "cols"),
          here::here("data", "ab_age.csv"))
saveRDS(files, list_file)

