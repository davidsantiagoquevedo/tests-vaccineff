library(arrow)
library(vaccineff)
library(ggplot2)
library(scales)
library(dplyr)
match_info <- function(data,
                       info_to_match) {
  matched_info <- unlist(
    tapply(data[[info_to_match]],
      data$subclass,
      function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          return(as.character(min(x, na.rm = TRUE)))
        }
      }
    )
  )
  # return data matched by subclass
  return(matched_info[data$subclass])
}

download <- FALSE

if (download) {
  url <- "https://github.com/mcastrolab/coronavac-ve-brazil/raw/main/data/persons_fortaleza_coronavac_cohort.parquet"
  file_name <- "data.parquet"
  dir.create("data/")
  file_path <- "data/"
  download.file(url, paste(file_path, file_name, sep = ""), mode = "wb")
}

df_raw <- read_parquet("data/data.parquet")
print(names(df_raw))

# Covid Deaths
df_raw$covid_death_date <- as.Date(df_raw$covid_death_date)
df_raw$death_date_not_covid <- as.Date(df_raw$death_date_not_covid)
df_raw$dose_two_date <- as.Date(df_raw$dose_two_date)

# Dates of study
start_cohort <- as.Date("2021-01-01")
end_cohort <- as.Date("2021-12-31")

# Truncate data_set

## deaths out of study interval
df_trun <- df_raw[(df_raw$covid_death_date <= end_cohort) |
  is.na(df_raw$covid_death_date),
]
df_trun <- df_trun[(df_trun$covid_death_date >= start_cohort) |
  is.na(df_trun$covid_death_date),
]

nrow(df_trun)
nrow(df_raw)

## vaccines out of study interval
mx0 <- max(df_trun$dose_two_date, na.rm = TRUE)
mi0 <- min(df_trun$dose_two_date, na.rm = TRUE)

df_trun$dose_two_date_trun <- as.Date(
  ifelse(df_trun$dose_two_date <= end_cohort,
    as.character(df_trun$dose_two_date),
    as.Date(NA_character_)
  )
)

mxf <- max(df_trun$dose_two_date_trun, na.rm = TRUE)
mif <- min(df_trun$dose_two_date_trun, na.rm = TRUE)

# Immunization date
df_trun$immunization <-
  get_immunization_date(
    data = df_trun,
    outcome_date_col = "covid_death_date",
    outcome_delay = 0,
    immunization_delay = 14,
    vacc_date_col = c("dose_two_date_trun"),
    end_cohort = end_cohort,
    take_first = FALSE
  )

df_trun$vaccine_status <- set_status(
  data = df_trun,
  col_names = "immunization",
  status = c("v", "u")
)

plot_coverage(
  data = df_trun,
  vacc_date_col = "immunization",
  unit = "month"
)

# Matching
# sample 100000
df_sample <- sample_n(df_trun, 100000)
matched_cohort <- match_cohort(data = df_sample,
  status_vacc_col = "vaccine_status",
  nearest = c(age = 1)
)

matched_cohort <- as.data.frame(matched_cohort)

head(matched_cohort[order(matched_cohort$subclass), ], n = 4)

censoring_date <- "death_date_not_covid"
matched_cohort$censoring_date <-  get_censoring_date_match(
  data = matched_cohort,
  outcome_date_col = "covid_death_date",
  censoring_date_col = censoring_date
)

matched_cohort$t0_fu <- as.Date(match_info(data = matched_cohort,
  info_to_match = "immunization"
))

matched_cohort$time_to_death <- get_time_to_event(
  data = matched_cohort,
  outcome_date_col = "covid_death_date",
  censoring_date_col = "censoring_date",
  start_cohort = start_cohort,
  end_cohort = end_cohort,
  start_from_immunization = TRUE,
  immunization_date_col = "t0_fu"
)

matched_cohort$matched_time <- match_info(data = matched_cohort,
  info_to_match = "time_to_death"
)
neg_times <- matched_cohort[matched_cohort$matched_time < 0, ]
head(neg_times[order(neg_times$subclass), ], n = 4)


cohort_corrected <- matched_cohort[matched_cohort$matched_time > 0, ]
print(paste0("Number of matched registers: ", nrow(matched_cohort)))
print(paste0("Number of registers with correct time-to-event: ", nrow(cohort_corrected)))

cohort_corrected$death_status <- set_status(
  data = cohort_corrected,
  col_names = "covid_death_date"
)

cohort_corrected$death_status <- ifelse(
  !is.na(cohort_corrected$censoring_date) &
    !is.na(cohort_corrected$covid_death_date) &
    (cohort_corrected$censoring_date <= cohort_corrected$covid_death_date),
  yes = 0,
  no = cohort_corrected$death_status
)


cohort_corrected$vaccine_status <- as.character(cohort_corrected$vaccine_status)

plt <- plot_survival(
  data = cohort_corrected,
  outcome_status_col = "death_status",
  time_to_event_col = "time_to_death",
  vacc_status_col = "vaccine_status",
  vaccinated_status = "v",
  unvaccinated_status = "u",
  vaccinated_color = "steelblue",
  unvaccinated_color = "darkred",
  start_cohort = start_cohort,
  end_cohort = end_cohort,
  percentage = TRUE,
  cumulative = TRUE
)
plt

surv <- plt$data
surv[surv$time == max(surv$time), ]
surv$loglog <- log(-log(surv$surv))

surv$logtime <- log(surv$time)
plt_loglog <- ggplot(data = surv) +
geom_step(ggplot2::aes(x = .data$logtime,
                       y = .data$loglog,
                       color = .data$strata)
  ) +
  theme_classic() +
  labs(x = "Log[Time to event] (Days)",
       y = "Log[-Log[Survival probability]]") +
  labs(colour = "Vaccine Status") +
  scale_color_manual(
    name = "Vaccine Status",
    values = c("steelblue", "darkred"),
    labels = c("v", "u")
  )
plt_loglog


p_thr <- 0.05

cx <- survival::coxph(  # nolint
  survival::Surv(time_to_death, death_status) ~
    vaccine_status , data = cohort_corrected
)

# Test the Proportional Hazards Assumption
test <- survival::cox.zph(cx)
hr <- round(exp(stats::coef(cx)), digits = 4)

# extract first and second element as limits
ci025 <- round(exp(stats::confint(cx)), 4)[1]
ci975 <- round(exp(stats::confint(cx)), 4)[2]

p <- test$table["GLOBAL", "p"]

df_summ <- data.frame(
  HR = hr,
  HR_low = ci025,
  HR_high = ci975,
  V_eff = 1 - hr,
  V_eff_low = 1 - ci975,
  V_eff_high = 1 - ci025,
  p_value = p # p_value must be a numeric
)

print(df_summ)