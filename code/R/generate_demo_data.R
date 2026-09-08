options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(tidyr)
})

SCRIPT_ROOT <- Sys.getenv("SEPSIS_PACO2_SCRIPT_ROOT", unset = "")
if (!nzchar(SCRIPT_ROOT)) stop("Set SEPSIS_PACO2_SCRIPT_ROOT before generating demo data.", call. = FALSE)
source(file.path(SCRIPT_ROOT, "00_functions", "release_config.R"))
paths <- release_paths()
if (!identical(paths$mode, "demo")) {
  stop("Synthetic data generation is available only in demo mode.", call. = FALSE)
}

set.seed(paths$seed)
cohorts <- c("MIMIC", "AmsterdamUMCdb", "Chinese cohort")
stays_per_cohort <- 500L
stay_index <- seq_len(stays_per_cohort * length(cohorts))
analysis_cohort <- rep(cohorts, each = stays_per_cohort)
cohort_shift <- unname(c(MIMIC = -0.5, AmsterdamUMCdb = 0.5, `Chinese cohort` = 1.0)[analysis_cohort])

baseline <- tibble(
  global_stay_id = sprintf("demo_%06d", stay_index),
  global_patient_id = sprintf("demo_patient_%06d", stay_index),
  global_hosp_id = sprintf("demo_hospital_%06d", stay_index),
  source_stay_id = sprintf("synthetic_stay_%06d", stay_index),
  source_patient_id = sprintf("synthetic_patient_%06d", stay_index),
  source_hosp_id = sprintf("synthetic_hospital_%06d", stay_index),
  source_database = paste0("synthetic_", tolower(gsub("[^A-Za-z]", "", analysis_cohort))),
  analysis_cohort = analysis_cohort,
  included_final = TRUE,
  adult_flag = TRUE,
  first_hosp_flag = TRUE,
  first_icu_flag = TRUE,
  icu_los_ge_24h = TRUE,
  sepsis3_window_flag = TRUE,
  outcome_28d_available = TRUE,
  missing_observable_day_paco2_flag = FALSE,
  age = pmin(90, pmax(18, round(rnorm(length(stay_index), 62 + cohort_shift, 14), 1))),
  sex = if_else(runif(length(stay_index)) < 0.55, "male", "female"),
  bmi = pmin(45, pmax(16, round(rnorm(length(stay_index), 26, 4.5), 1))),
  cci = pmax(0L, rpois(length(stay_index), 2)),
  copd = runif(length(stay_index)) < 0.15,
  baseline_sofa = pmin(22, pmax(2, round(rnorm(length(stay_index), 10, 3), 1))),
  baseline_lactate = pmin(15, pmax(0.4, round(rlnorm(length(stay_index), log(2.5), 0.55), 2))),
  baseline_ph = pmin(7.60, pmax(6.90, round(rnorm(length(stay_index), 7.34, 0.08), 3))),
  baseline_hco3 = pmin(40, pmax(6, round(rnorm(length(stay_index), 22, 5), 1))),
  baseline_pao2 = pmin(450, pmax(35, round(rlnorm(length(stay_index), log(115), 0.45), 1))),
  baseline_fio2 = pmin(1, pmax(0.21, round(runif(length(stay_index), 0.28, 0.85), 2))),
  baseline_mv = runif(length(stay_index)) < 0.70,
  baseline_vasopressor = runif(length(stay_index)) < 0.55,
  baseline_crrt_rrt = runif(length(stay_index)) < 0.10,
  icu_los_days = round(runif(length(stay_index), 7, 14), 2)
) %>%
  mutate(
    baseline_pf_ratio = round(.data$baseline_pao2 / .data$baseline_fio2, 1),
    mean_demo_paco2 = pmin(65, pmax(25, rnorm(n(), 41 + cohort_shift, 7))),
    demo_risk = plogis(-2.0 + 0.12 * pmax(0, 35 - .data$mean_demo_paco2) + 0.10 * pmax(0, .data$mean_demo_paco2 - 50) + 0.05 * (.data$baseline_sofa - 10)),
    death_28d = runif(n()) < .data$demo_risk,
    time_to_event_28d_days = if_else(.data$death_28d, runif(n(), 1.25, 27.75), 28),
    death_time = .data$time_to_event_28d_days
  )

day_long <- crossing(
  global_stay_id = baseline$global_stay_id,
  icu_day = seq_len(7L)
) %>%
  left_join(baseline, by = "global_stay_id") %>%
  mutate(
    source_icu_day = .data$icu_day,
    observable_icu_day = TRUE,
    landmark_eligible = .data$time_to_event_28d_days > .data$icu_day,
    time_from_window_end_to_event_or_censor_days = pmax(.data$time_to_event_28d_days - .data$icu_day, 0),
    n_paco2 = sample(1:4, n(), replace = TRUE),
    twa_paco2 = pmin(70, pmax(20, round(.data$mean_demo_paco2 + 0.35 * (.data$icu_day - 1) + rnorm(n(), 0, 2.5), 2))),
    daily_sofa = pmin(24, pmax(0, round(.data$baseline_sofa + rnorm(n(), -0.15 * (.data$icu_day - 1), 1.6), 1))),
    daily_lactate = pmin(18, pmax(0.3, round(.data$baseline_lactate * exp(rnorm(n(), -0.04 * (.data$icu_day - 1), 0.25)), 2))),
    daily_mv = .data$baseline_mv | runif(n()) < 0.08,
    daily_vasopressor = .data$baseline_vasopressor & runif(n()) < pmax(0.15, 0.92 - 0.08 * (.data$icu_day - 1)),
    daily_crrt_rrt = .data$baseline_crrt_rrt | runif(n()) < 0.03,
    daily_ph = pmin(7.60, pmax(6.90, round(.data$baseline_ph - 0.003 * (.data$twa_paco2 - 41) + rnorm(n(), 0, 0.035), 3))),
    daily_hco3 = pmin(42, pmax(5, round(.data$baseline_hco3 + rnorm(n(), 0, 2.3), 1))),
    daily_pao2 = pmin(500, pmax(30, round(.data$baseline_pao2 * exp(rnorm(n(), 0, 0.18)), 1))),
    daily_fio2 = pmin(1, pmax(0.21, round(.data$baseline_fio2 + rnorm(n(), 0, 0.06), 2))),
    daily_pf_ratio = round(.data$daily_pao2 / .data$daily_fio2, 1),
    daily_oxygenation_value = .data$daily_pf_ratio,
    daily_oxygenation_type = "pf_ratio"
  ) %>%
  select(
    source_database, analysis_cohort, source_patient_id, source_hosp_id,
    source_stay_id, global_patient_id, global_hosp_id, global_stay_id,
    icu_day, source_icu_day, observable_icu_day, landmark_eligible,
    time_from_window_end_to_event_or_censor_days, n_paco2,
    twa_paco2, daily_sofa, daily_lactate, daily_mv, daily_vasopressor,
    daily_crrt_rrt, daily_ph, daily_hco3, daily_pao2, daily_fio2,
    daily_pf_ratio, daily_oxygenation_value, daily_oxygenation_type,
    included_final, age, sex, bmi, cci, copd, baseline_sofa, baseline_lactate,
    baseline_ph, baseline_hco3, baseline_pao2, baseline_fio2, baseline_pf_ratio,
    baseline_mv, baseline_vasopressor, baseline_crrt_rrt, death_28d,
    time_to_event_28d_days, death_time, icu_los_days
  )

raw_paco2 <- day_long %>%
  transmute(
    global_stay_id, analysis_cohort, icu_day, paco2 = .data$twa_paco2,
    used_for_twa = TRUE
  )

pooled_root <- file.path(paths$input_root, "pooled")
dir.create(pooled_root, recursive = TRUE, showWarnings = FALSE)
write_parquet(baseline %>% select(-mean_demo_paco2, -demo_risk), file.path(pooled_root, "all_baseline_outcome.parquet"))
write_parquet(day_long, file.path(pooled_root, "all_day_long.parquet"))
write_parquet(raw_paco2, file.path(pooled_root, "all_raw_paco2.parquet"))

write_csv(
  tibble(
    item = c("mode", "seed", "synthetic_stays", "synthetic_stay_days", "synthetic_deaths_28d"),
    value = c(paths$mode, paths$seed, nrow(baseline), nrow(day_long), sum(baseline$death_28d))
  ),
  file.path(paths$input_root, "demo_generation_metadata.csv")
)

message("Generated fully simulated demo data in: ", paths$input_root)
