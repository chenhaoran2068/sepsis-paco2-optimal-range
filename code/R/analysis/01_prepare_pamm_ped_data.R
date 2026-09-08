options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(readr)
})

SCRIPT_ROOT <- Sys.getenv("SEPSIS_PACO2_SCRIPT_ROOT", unset = "")
if (!nzchar(SCRIPT_ROOT)) stop("Set SEPSIS_PACO2_SCRIPT_ROOT before running this release-candidate script.", call. = FALSE)
source(file.path(SCRIPT_ROOT, "00_functions", "release_config.R"))
paths <- release_paths()
PROJECT_ROOT <- paths$release_root
RERUN_ID <- paths$run_id
ANALYSIS_DATA_ROOT <- paths$input_root
RESULT_ROOT <- paths$result_root
SCRIPT_ROOT <- paths$script_root
ANALYSIS_NTHREADS <- paths$nthreads
set.seed(paths$seed)
POOLED_DIR <- file.path(ANALYSIS_DATA_ROOT, "pooled")
OUT_DIR <- paths$pamm_root
QC_DIR <- file.path(OUT_DIR, "_qc")
MAX_EXPOSURE_DAY <- 7L
MAX_FOLLOWUP_DAY <- 28L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

message("Preparing PAMM/PED data...")

write_qc <- function(data, file_name) {
  write_csv(data, file.path(QC_DIR, file_name), na = "")
}

baseline <- read_parquet(file.path(POOLED_DIR, "all_baseline_outcome.parquet")) %>%
  filter(.data$included_final) %>%
  mutate(
    sex = case_when(
      tolower(.data$sex) %in% c("male", "m") ~ "male",
      tolower(.data$sex) %in% c("female", "f") ~ "female",
      TRUE ~ NA_character_
    ),
    followup_stop_days = if_else(
      .data$death_28d,
      pmin(.data$time_to_event_28d_days, MAX_FOLLOWUP_DAY),
      MAX_FOLLOWUP_DAY
    ),
    followup_stop_days = pmax(0, pmin(MAX_FOLLOWUP_DAY, .data$followup_stop_days)),
    n_ped_intervals = pmax(1L, as.integer(ceiling(.data$followup_stop_days)))
  )

day_long <- read_parquet(file.path(POOLED_DIR, "all_day_long.parquet")) %>%
  filter(.data$included_final, .data$observable_icu_day) %>%
  select(
    source_database,
    analysis_cohort,
    global_patient_id,
    global_hosp_id,
    global_stay_id,
    source_icu_day = icu_day,
    twa_paco2,
    n_paco2,
    daily_sofa,
    daily_lactate,
    daily_mv,
    daily_vasopressor,
    daily_crrt_rrt,
    daily_ph,
    daily_hco3,
    daily_pao2,
    daily_fio2,
    daily_pf_ratio,
    daily_oxygenation_value,
    daily_oxygenation_type
  )

last_observed_day <- day_long %>%
  group_by(.data$global_stay_id) %>%
  summarise(
    first_observable_icu_day = min(.data$source_icu_day, na.rm = TRUE),
    last_observable_icu_day = max(.data$source_icu_day, na.rm = TRUE),
    n_observable_icu_days = n_distinct(.data$source_icu_day),
    .groups = "drop"
  )

missing_observable <- anti_join(
  baseline %>% select(.data$global_stay_id),
  last_observed_day,
  by = "global_stay_id"
)

if (nrow(missing_observable) > 0) {
  stop("Some final included stays have no observable ICU day rows. See upstream pooled QC.", call. = FALSE)
}

ped_intervals <- baseline %>%
  left_join(last_observed_day, by = "global_stay_id") %>%
  uncount(.data$n_ped_intervals, .id = "ped_interval") %>%
  mutate(
    tstart = as.numeric(.data$ped_interval - 1L),
    tend = pmin(as.numeric(.data$ped_interval), .data$followup_stop_days),
    interval_length_days = pmax(.data$tend - .data$tstart, 0),
    ped_status = .data$death_28d &
      .data$time_to_event_28d_days > .data$tstart &
      .data$time_to_event_28d_days <= .data$tend,
    offset_log_interval = log(pmax(.data$interval_length_days, 1e-08)),
    followup_day_mid = (.data$tstart + .data$tend) / 2,
    nominal_followup_day = pmin(
      MAX_FOLLOWUP_DAY,
      pmax(1L, as.integer(floor(.data$tstart) + 1L))
    ),
    nominal_exposure_day = pmin(MAX_EXPOSURE_DAY, .data$nominal_followup_day),
    source_icu_day = pmin(.data$nominal_exposure_day, .data$last_observable_icu_day),
    covariate_carried_forward = .data$source_icu_day != .data$nominal_followup_day,
    post_last_observable_icu_day = .data$nominal_followup_day > .data$last_observable_icu_day,
    post_exposure_window_day7 = .data$nominal_followup_day > MAX_EXPOSURE_DAY,
    in_icu_at_interval_start = .data$tstart < .data$icu_los_days,
    after_icu_discharge_interval = !.data$in_icu_at_interval_start,
    ped_time_scale = "days_after_icu_admission",
    exposure_assignment_rule = if_else(
      .data$covariate_carried_forward,
      "last_observable_icu_day_carried_forward",
      "same_followup_day"
    )
  ) %>%
  filter(.data$interval_length_days > 0)

ped_full <- ped_intervals %>%
  left_join(
    day_long %>% select(-source_database, -analysis_cohort, -global_patient_id, -global_hosp_id),
    by = c("global_stay_id", "source_icu_day")
  ) %>%
  mutate(
    daily_mv = coalesce(.data$daily_mv, FALSE),
    daily_vasopressor = coalesce(.data$daily_vasopressor, FALSE),
    daily_crrt_rrt = coalesce(.data$daily_crrt_rrt, FALSE),
    twa_paco2_per10 = .data$twa_paco2 / 10,
    twa_paco2_centered_40 = .data$twa_paco2 - 40,
    age_per10 = .data$age / 10,
    bmi_per5 = .data$bmi / 5,
    daily_lactate_log1p = log1p(.data$daily_lactate),
    daily_oxygenation_per100 = .data$daily_oxygenation_value / 100,
    complete_model1 = !is.na(.data$twa_paco2),
    complete_model2 = .data$complete_model1 &
      !is.na(.data$analysis_cohort) &
      !is.na(.data$age) &
      !is.na(.data$sex) &
      !is.na(.data$bmi),
    complete_model3 = .data$complete_model2 &
      !is.na(.data$daily_sofa) &
      !is.na(.data$daily_mv) &
      !is.na(.data$daily_vasopressor) &
      !is.na(.data$daily_crrt_rrt) &
      !is.na(.data$daily_lactate) &
      !is.na(.data$daily_hco3) &
      !is.na(.data$daily_oxygenation_value) &
      !is.na(.data$daily_oxygenation_type),
    complete_sensitivity_ph = .data$complete_model3 & !is.na(.data$daily_ph)
  ) %>%
  select(
    source_database,
    analysis_cohort,
    source_patient_id,
    source_hosp_id,
    source_stay_id,
    global_patient_id,
    global_hosp_id,
    global_stay_id,
    ped_interval,
    tstart,
    tend,
    interval_length_days,
    offset_log_interval,
    followup_day_mid,
    nominal_followup_day,
    nominal_exposure_day,
    source_icu_day,
    first_observable_icu_day,
    last_observable_icu_day,
    n_observable_icu_days,
    covariate_carried_forward,
    post_last_observable_icu_day,
    post_exposure_window_day7,
    in_icu_at_interval_start,
    after_icu_discharge_interval,
    ped_time_scale,
    exposure_assignment_rule,
    ped_status,
    death_28d,
    death_time,
    time_to_event_28d_days,
    followup_stop_days,
    age,
    age_per10,
    sex,
    bmi,
    bmi_per5,
    cci,
    copd,
    baseline_sofa,
    baseline_lactate,
    baseline_mv,
    baseline_vasopressor,
    baseline_crrt_rrt,
    twa_paco2,
    twa_paco2_per10,
    twa_paco2_centered_40,
    n_paco2,
    daily_sofa,
    daily_lactate,
    daily_lactate_log1p,
    daily_mv,
    daily_vasopressor,
    daily_crrt_rrt,
    daily_hco3,
    daily_oxygenation_value,
    daily_oxygenation_per100,
    daily_oxygenation_type,
    daily_ph,
    daily_pao2,
    daily_fio2,
    daily_pf_ratio,
    complete_model1,
    complete_model2,
    complete_model3,
    complete_sensitivity_ph
  )

model_variable_sets <- tibble(
  model_id = c(
    "model1",
    "model2",
    "model3",
    "sensitivity_ph",
    "retained_not_primary_pooled"
  ),
  role = c(
    "unadjusted PAMM input",
    "time-fixed adjusted PAMM input",
    "primary adjusted PAMM input",
    "pH-adjusted sensitivity PAMM input",
    "available for cohort-specific or secondary checks"
  ),
  variables = c(
    "tstart;tend;ped_status;offset_log_interval;followup_day_mid;twa_paco2",
    "model1 variables;analysis_cohort;age;sex;bmi",
    "model2 variables;daily_sofa;daily_mv;daily_vasopressor;daily_crrt_rrt;daily_lactate;daily_hco3;daily_oxygenation_value;daily_oxygenation_type",
    "model3 variables;daily_ph",
    "cci;copd;baseline_sofa;baseline_lactate;baseline_mv;baseline_vasopressor;baseline_crrt_rrt;daily_pao2;daily_fio2;daily_pf_ratio"
  ),
  complete_case_flag = c(
    "complete_model1",
    "complete_model2",
    "complete_model3",
    "complete_sensitivity_ph",
    NA_character_
  ),
  note = c(
    "Exposure-only model; baseline hazard and nonlinear PaCO2 terms are specified at fitting stage.",
    "CCI and COPD are retained but not included in the current pooled Model 2; CCI remains unavailable for Amsterdam/Xuzhou, while COPD is reserved for baseline reporting and subgroup work unless the user approves primary adjustment.",
    "First PED version uses complete-case rows for daily covariates; multiple imputation is not performed at this step.",
    "pH is added only here because it may be on the causal pathway from PaCO2 to mortality.",
    "These fields remain in the PED table for audit, Table S1, cohort-specific checks, subgroup work, or later user-approved revisions."
  )
)

model_input_rules <- tibble(
  item = c(
    "source_tables",
    "analysis_population",
    "followup_time_scale",
    "ped_interval_length",
    "outcome",
    "exposure",
    "time_varying_covariate_source",
    "post_icu_or_post_day7_assignment",
    "primary_complete_case_policy",
    "pH_policy",
    "comorbidity_policy",
    "apache_policy"
  ),
  value = c(
    "data/data_for_analysis/pooled/all_baseline_outcome.parquet and all_day_long.parquet",
    "included_final == TRUE",
    "days after ICU admission, censored at 28 days",
    "one row per patient-day interval; final interval truncated at death time if death occurs before 28 days",
    "ped_status indicates death within the current interval; death_28d remains the patient-level 28-day outcome",
    "daily TWA-PaCO2 from ICU Day 1-7 observable windows",
    "same ICU day while available",
    "after ICU Day 7 or after the last observable ICU day, the last observable ICU-day covariate set is carried forward and flagged",
    "first PED version writes model-specific complete-case datasets; imputation is not performed here",
    "daily pH is excluded from Model 3 and included only in ped_sensitivity_ph",
    "CCI is retained for audit but not used in pooled Model 2/3 because it remains unavailable for Amsterdam and Xuzhou; COPD is retained for baseline reporting and subgroup work, not current primary adjustment",
    "APACHE/APS/SAPS/OASIS are not used"
  )
)

write_parquet(ped_full, file.path(OUT_DIR, "ped_full_28d_locf.parquet"))
write_parquet(ped_full %>% filter(.data$complete_model1), file.path(OUT_DIR, "ped_model1.parquet"))
write_parquet(ped_full %>% filter(.data$complete_model2), file.path(OUT_DIR, "ped_model2.parquet"))
write_parquet(ped_full %>% filter(.data$complete_model3), file.path(OUT_DIR, "ped_model3.parquet"))
write_parquet(ped_full %>% filter(.data$complete_sensitivity_ph), file.path(OUT_DIR, "ped_sensitivity_ph.parquet"))

ped_qc_specs <- tibble(
  model_id = c("ped_full_28d_locf", "ped_model1", "ped_model2", "ped_model3", "ped_sensitivity_ph"),
  complete_case_flag = c(NA_character_, "complete_model1", "complete_model2", "complete_model3", "complete_sensitivity_ph")
)

ped_qc_by_model <- bind_rows(lapply(seq_len(nrow(ped_qc_specs)), function(i) {
  flag <- ped_qc_specs$complete_case_flag[i]
  model_data <- if (is.na(flag)) {
    ped_full
  } else {
    ped_full %>% filter(.data[[flag]])
  }

  tibble(
    model_id = ped_qc_specs$model_id[i],
    complete_case_flag = flag,
    n_rows = nrow(model_data),
    n_unique_stays = n_distinct(model_data$global_stay_id),
    n_events = sum(model_data$ped_status, na.rm = TRUE)
  )
}))

ped_qc_by_model_cohort <- bind_rows(
  ped_full %>% mutate(model_id = "ped_full_28d_locf", in_model = TRUE),
  ped_full %>% mutate(model_id = "ped_model1", in_model = .data$complete_model1),
  ped_full %>% mutate(model_id = "ped_model2", in_model = .data$complete_model2),
  ped_full %>% mutate(model_id = "ped_model3", in_model = .data$complete_model3),
  ped_full %>% mutate(model_id = "ped_sensitivity_ph", in_model = .data$complete_sensitivity_ph)
) %>%
  filter(.data$in_model) %>%
  group_by(.data$model_id, .data$analysis_cohort) %>%
  summarise(
    n_rows = n(),
    n_unique_stays = n_distinct(.data$global_stay_id),
    n_events = sum(.data$ped_status, na.rm = TRUE),
    n_carried_forward_rows = sum(.data$covariate_carried_forward, na.rm = TRUE),
    n_post_day7_rows = sum(.data$post_exposure_window_day7, na.rm = TRUE),
    .groups = "drop"
  )

ped_missingness_full <- ped_full %>%
  summarise(
    n_rows = n(),
    across(
      c(
        twa_paco2,
        age,
        sex,
        bmi,
        cci,
        copd,
        daily_sofa,
        daily_lactate,
        daily_hco3,
        daily_oxygenation_value,
        daily_oxygenation_type,
        daily_ph
      ),
      ~ sum(is.na(.x)),
      .names = "missing_{.col}"
    ),
    .by = .data$analysis_cohort
  )

ped_integrity_checks <- tibble(
  check_name = c(
    "baseline_final_patients",
    "ped_full_unique_stays",
    "ped_full_rows_positive_interval",
    "ped_full_duplicate_stay_interval",
    "ped_status_total_events",
    "baseline_death_28d_total",
    "ped_event_matches_baseline_death_28d",
    "model1_missing_twa_rows",
    "model3_includes_ph",
    "sensitivity_ph_missing_ph_rows",
    "source_icu_day_outside_1_to_7",
    "zero_interval_rows"
  ),
  value = c(
    nrow(baseline),
    n_distinct(ped_full$global_stay_id),
    nrow(ped_full),
    anyDuplicated(paste(ped_full$global_stay_id, ped_full$ped_interval)),
    sum(ped_full$ped_status, na.rm = TRUE),
    sum(baseline$death_28d, na.rm = TRUE),
    as.integer(sum(ped_full$ped_status, na.rm = TRUE) == sum(baseline$death_28d, na.rm = TRUE)),
    sum(ped_full$complete_model1 & is.na(ped_full$twa_paco2), na.rm = TRUE),
    sum(ped_full$complete_model3 & !is.na(ped_full$daily_ph), na.rm = TRUE),
    sum(ped_full$complete_sensitivity_ph & is.na(ped_full$daily_ph), na.rm = TRUE),
    sum(is.na(ped_full$source_icu_day) | ped_full$source_icu_day < 1L | ped_full$source_icu_day > MAX_EXPOSURE_DAY, na.rm = TRUE),
    sum(ped_full$interval_length_days <= 0, na.rm = TRUE)
  )
)

write_qc(model_variable_sets, "model_variable_sets.csv")
write_qc(model_input_rules, "model_input_rules.csv")
write_qc(ped_qc_by_model, "ped_qc_by_model.csv")
write_qc(ped_qc_by_model_cohort, "ped_qc_by_model_cohort.csv")
write_qc(ped_missingness_full, "ped_missingness_full.csv")
write_qc(ped_integrity_checks, "ped_integrity_checks.csv")

metadata <- tibble(
  item = c(
    "output_dir",
    "build_time",
    "ped_time_scale",
    "max_followup_day",
    "max_exposure_day",
    "primary_ped_file",
    "primary_adjusted_file",
    "sensitivity_ph_file"
  ),
  value = c(
    OUT_DIR,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    "days after ICU admission",
    as.character(MAX_FOLLOWUP_DAY),
    as.character(MAX_EXPOSURE_DAY),
    file.path(OUT_DIR, "ped_full_28d_locf.parquet"),
    file.path(OUT_DIR, "ped_model3.parquet"),
    file.path(OUT_DIR, "ped_sensitivity_ph.parquet")
  )
)

write_qc(metadata, "build_metadata.csv")

message("Done.")
message("Output directory: ", OUT_DIR)
message("PED full rows: ", nrow(ped_full))
message("PED Model 3 rows: ", sum(ped_full$complete_model3, na.rm = TRUE))
