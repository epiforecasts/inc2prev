# Code adapted from: https://github.com/cmmid/pcr-profile/blob/main/run_analysis.R
# Adaptions are code cleaning and tidying up of the merged data.frame
# Paper: https://doi.org/10.1186/s12916-021-01982-x

library(data.table)
library(ggplot2)
library(ggnewscale)

test_final <- data.table::fread("https://raw.githubusercontent.com/cmmid/pcr-profile/main/test_data.csv") # nolint
symp_final <- data.table::fread("https://raw.githubusercontent.com/cmmid/pcr-profile/main/symptom_data.csv") # nolint
symp_final$date <- as.Date(symp_final$date)

## Set days relative to a chosen start date
start_date <- as.Date(min(test_final$date) - 10)
test_final[, date := as.Date(date)]
test_final[, serology_date := as.Date(serology_date)]
test_final[, day := as.integer(date - start_date)]
test_final[, serology_day := as.integer(serology_date - start_date)]

## Add initial asymptomatic reports on enrollment day
symp_final <- rbind(
  symp_final, test_final[, .(date = min(date), symptom = FALSE), by = num_id]
)

## Find first symptomatic report dates and last asymptomatic report dates
symp_final[,
  first_symp := min(date[symptom == TRUE]), by = num_id
]
symp_final[,
 last_asym := max(date[symptom == FALSE & date < first_symp]), by = num_id
]

first_last_df <- symp_final[, .(
  first_sym_day = unique(as.integer(first_symp - start_date)),
  last_asym_day = unique(as.integer(last_asym - start_date))),
  by = num_id
]

symp_final[, day := as.integer(date - start_date)][, date := NULL]

setkey(first_last_df, num_id)
setkey(test_final, num_id)

test_final <- merge(test_final, first_last_df)

# Merge testing and symptom data into one data table
pcr_testing <- merge(
  symp_final, test_final, by = c("num_id", "day"), all.y = TRUE
)

# Add infection day upper bound
# Must occur prior to symptoms or first positive  test
pcr_testing[,
  inf_upper_bound := min(
    day[pcr_result == TRUE], first_sym_day, na.rm = TRUE
    ),
  by = num_id
]

# Save tidy data
fwrite(
  pcr_testing,
  here("data", "pcr_testing.csv")
)


# Check data looks like the paper using figure 1 code
figure1 <- function(dfy = NULL) {

  cols1 <- scales::viridis_pal(option = "D")(10)
  cols2 <- scales::viridis_pal(option = "A")(10)

  # Generate figure 1
  fig1 <- dfy |>
    ggplot(aes(x = day  - last_asym_day, y = as.factor(num_id))) +
    geom_point(
      aes(fill = symptom), size = 3, shape = 21, stroke = 1.5, col = "white"
    ) +
    geom_point(aes(col = pcr_result), size = 3, shape = 21, stroke = 1.5) +
    theme_bw() +
    scale_color_manual(
      values = cols1[c(2, 7)], name = "PCR result",
      labels = c("Negative", "Positive")
    ) +
    scale_fill_manual(
      values = c("white","red"), name = "Symptoms", na.value = "grey40"
    ) +
    labs(x = "Days since last asymptomatic report", y = "Participant ID") +
    scale_x_continuous(
      minor_breaks = seq(-19, 38, 1), breaks = seq(-16, 36, 4)
    ) +
    theme(panel.grid.minor = element_line(size = (0.2), colour = "grey")) +
    coord_cartesian(xlim = c(-19, 38)) +
    new_scale_color() +
    geom_point(
      data = dfy[!is.na(serology_day),
        .(serology_day = serology_day[1] - last_asym_day), num_id][,
        sero := TRUE],
        inherit.aes = FALSE,
        aes(x = serology_day, y = as.factor(num_id), col = sero),
        pch = 4, size = 4
    ) +
    scale_color_manual(
      name = "Serology result", values = "black", labels = "Negative"
  )
  return(fig1)
}

fig1 <- figure1(pcr_testing)

# Save Figure 1
ggsave(
  fig1, filename = "data-raw/figure1.png",
  height = 20, width = 40, units = "cm"
)


pcr_test_data_to_list <- function(dt) {

if (!is.null(dt)) {
  ids <- unique(dt[,
   .(num_id, first_sym_day, last_asym_day, inf_upper_bound)])

  dat <- list(
    pcr_n = nrow(ids),
    pcr_p = nrow(dt),
    pcr_id = ids$num_id,
    pcr_test_day = dt$day,
    pcr_result = as.numeric(dt$pcr_result),
    pcr_sym_at_test = ids$first_sym_day,
    pcr_last_asym_at_test = ids$last_asym_day,
    pcr_inf_upper_bound = ids$inf_upper_bound
  )
}else{
  dat <- list(
    pcr_n = 0,
    pcr_p = 1,
    pcr_test_day = 1,
    pcr_result = 1,
    pcr_sym_at_test = 1,
    pcr_last_asym_at_test = 1,
    pcr_inf_upper_bound = 1
  )
}
  return(dat)
}
