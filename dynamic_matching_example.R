library("vaccineff")
data(cohortdata)

new_calendar_match <- function(data,
                               matched_cohort,
                               rolling_date,
                               outcome_date,
                               censoring_date,
                               calendar_time,
                               exact = NULL,
                               nearest = NULL) {
  # Select new vaccinated people elegible for matching at rolling_date
  new_vaccinated <- data[
    (data[[rolling_date]] == calendar_time) &
      !is.na(data[[rolling_date]]),
  ]
  # Set vaccinated status
  new_vaccinated$vacc_status <- "v"
  # Select unvaccinated people elegible for matching at rolling_date
  # 1. vaccine date must be greater than calendar time or NA
  unvaccinated <- data[
    (data[[rolling_date]] > calendar_time) |
      is.na(data[[rolling_date]]),
  ]
  print(nrow(unvaccinated))
  # 2. Not previously matched
  unvaccinated <- unvaccinated[
    !(unvaccinated$id %in% matched_cohort$id),
  ]
  print(nrow(unvaccinated))
  # 3. Outcome date greater than calendar time or NA
  unvaccinated <- unvaccinated[
    (unvaccinated[[outcome_date]] > calendar_time) |
      is.na(unvaccinated[[outcome_date]]),
  ]
  print(nrow(unvaccinated))
  # 4. censoring date greater than calendar time or NA
  unvaccinated <- unvaccinated[
    (unvaccinated[[censoring_date]] > calendar_time) |
      is.na(unvaccinated[[censoring_date]]),
  ]
  print(nrow(unvaccinated))
  # 5. Set vaccinated status as unvaccinated by default
  unvaccinated$vacc_status <- "u"
  # Combine both data frames
  unmatched <- rbind(new_vaccinated, unvaccinated)
  new_matches <- match_cohort(data = unmatched,
    status_vacc_col = "vacc_status",
    exact = exact,
    nearest = nearest
  )
  # t0 follow up
  new_matches$t0_follow_up <- calendar_time
  # censoring dates for new matches
  new_matches$censoring_date_match <- get_censoring_date_match(
    data = new_matches,
    outcome_date_col = outcome_date,
    censoring_date_col = censoring_date
  )

  return(new_matches)
}

get_censoring_date_par_vacc <- function(matched_cohort,
                                        new_matches,
                                        outcome_date,
                                        censoring_date_match,
                                        calendar_time) {
  # censoring_date_par_vacc does not come in the first iteration of the
  # matching. It must be defined filled by NAs
  if (!("censoring_date_par_vacc" %in% names(matched_cohort))) {
    matched_cohort$censoring_date_par_vacc <- as.Date(NA)
  }
  # get subclasses of individuals vaccinated censored
  # when unvaccinated partner receives vaccine
  subclass_new_censored <- matched_cohort[
    (matched_cohort$id %in% new_matches$id),
  ]$subclass

  # calendar_time is assigned to censoring_date_par_vacc if
  # 1. the individual was previously vaccinated
  # 2. censoring_date_match is greater than calendar_time
  # 3. outcome_date is greater than calendar_time
  matched_cohort$censoring_date_par_vacc <-
    as.Date(
      ifelse((matched_cohort$subclass %in% subclass_new_censored) &
          (matched_cohort$vacc_status == "v") &
          ((matched_cohort[[censoring_date_match]] > calendar_time) |
             is.na(matched_cohort[[censoring_date_match]])) &
          ((matched_cohort[[censoring_date_match]] >
              matched_cohort[[outcome_date]]) |
             is.na(matched_cohort[[outcome_date]])),
        yes =  as.character(calendar_time),
        no = as.character(matched_cohort$censoring_date_par_vacc)
      )
    )
  return(matched_cohort$censoring_date_par_vacc)
}

data <- cohortdata
outcome_date <- "death_date"
censoring_date <- "death_other_causes"
immunization_delay <- 14
end_cohort <- as.Date("2044-12-31")
vacc_date_col <- "vaccine_date_2"

data$immunization <-
  get_immunization_date(
    data = data,
    outcome_date_col = outcome_date,
    outcome_delay = 0,
    immunization_delay = immunization_delay,
    vacc_date_col = vacc_date_col,
    end_cohort = end_cohort,
    take_first = FALSE
  )

rolling_date <- "immunization"
nearest <- c(age = 1)
matched_cohort <- data.frame()

calendar_time <- as.Date("2044-04-01")
new_matches <- new_calendar_match(data = data,
  matched_cohort = matched_cohort,
  rolling_date = rolling_date,
  outcome_date = outcome_date,
  censoring_date = censoring_date,
  calendar_time = calendar_time,
  nearest = nearest
)

# MatchIt generates subclass IDs starting from 1.
# To avoid repeating, it is necessary to add the lenght of the current
# dataset with the matched population
new_matches$subclass <- as.numeric(new_matches$subclass) + nrow(matched_cohort)

# If is the first iteration in time matched_cohort is directly updated
# by the new dataset matched
matched_cohort <- rbind(matched_cohort, new_matches)

# New match
calendar_time <- as.Date(calendar_time + 1)

new_matches <- new_calendar_match(data = data,
  matched_cohort = matched_cohort,
  rolling_date = rolling_date,
  outcome_date = outcome_date,
  censoring_date = censoring_date,
  calendar_time = calendar_time,
  nearest = nearest
)

matched_cohort$censoring_date_par_vacc <-
  get_censoring_date_par_vacc(matched_cohort = matched_cohort,
    new_matches = new_matches,
    outcome_date = outcome_date,
    censoring_date_match = "censoring_date_match",
    calendar_time = calendar_time
  )


# previously matched registers
id_pre_matched <-  matched_cohort[
  (matched_cohort$id %in% new_matches$id),
]$id

# extract old subclasses of new censored registers
subclass_new_censored <- matched_cohort[
  !is.na(matched_cohort$censoring_date_par_vacc) &
    matched_cohort$censoring_date_par_vacc == calendar_time,
]$subclass

# extract id's of accepted transitions for u to v
# a transitio is accepted if ex-partner is new censored
id_accepted_u_to_v <- matched_cohort[
  matched_cohort$subclass %in% subclass_new_censored &
    matched_cohort$vacc_status == "u",
]$id

# extract id's of rejected transitions for u to v
id_rejected_u_to_v <- id_pre_matched[
  !(id_pre_matched %in% id_accepted_u_to_v)
]

# keep track of censored registers by new vaccinated
new_censored <- matched_cohort[
  matched_cohort$subclass %in% subclass_new_censored
]

# remove old unvaccinated register of accepted transition u to v
matched_cohort <- matched_cohort[
  !(matched_cohort$id %in% id_accepted_u_to_v)
]

# remove new vaccinated register of rejected transition u to v
new_matches <- new_matches[
  !(new_matches$id %in% id_rejected_u_to_v)
]

matched_cohort <- rbind(matched_cohort, new_matches)


# Coverage manual check

coh_coverage(data = cohortdata,
  vacc_date_col = rolling_date,
  unit = "day"
)


table(new_vaccinated$immunization, useNA = "ifany")
head(sort(unvaccinated[!is.na(unvaccinated$immunization), ]$immunization))
