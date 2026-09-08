options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(openxlsx)
  library(readr)
  library(survival)
  library(tibble)
  library(tidyr)
})

# 配置环境（载入包、载入数据） ------------------------------------------------

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
PAMM_DIR <- paths$pamm_root

FIG_MAIN_DIR <- file.path(RESULT_ROOT, "figures", "main")
FIG_SUPP_DIR <- file.path(RESULT_ROOT, "figures", "supplementary")
TABLE_SUPP_DIR <- file.path(RESULT_ROOT, "tables", "supplementary")
MANUSCRIPT_DIR <- file.path(RESULT_ROOT, "manuscript")
MODEL_OUT_DIR <- file.path(PAMM_DIR, "result4_models")

dir.create(FIG_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Running Result 4: subgroup analyses and pH-adjusted sensitivity...")

PRIMARY_K <- 7L

cohort_display <- c(
  MIMIC = "MIMIC-IV",
  AmsterdamUMCdb = "AmsterdamUMCdb",
  `Chinese cohort` = "Chinese"
)

scope_display <- c(
  pooled = "Pooled",
  mimic = "MIMIC-IV",
  amsterdam = "AmsterdamUMCdb",
  chinese = "Chinese"
)

scope_ids <- names(scope_display)
scope_to_cohort <- c(
  mimic = "MIMIC",
  amsterdam = "AmsterdamUMCdb",
  chinese = "Chinese cohort"
)

figure_colors <- c(
  "Primary pH-available Model 3" = "#555555",
  "pH-adjusted Model 3" = "#D95F0E",
  "Any cumulative high-risk exposure" = "#2C7FB8"
)

base_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 1)
    )
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  names(sort(table(x), decreasing = TRUE))[1]
}

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = digits))
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "P < .001",
    TRUE ~ paste0("P = ", sub("^0", "", sprintf("%.3f", x)))
  )
}

fmt_p_plain <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "< .001",
    TRUE ~ sub("^0", "", sprintf("%.3f", x))
  )
}

fmt_estimate <- function(hr, low, high) {
  paste0(sprintf("%.2f", hr), "（95% CI ", sprintf("%.2f", low), "-", sprintf("%.2f", high), "）")
}

fmt_hr_ci_english <- function(hr, low, high) {
  ifelse(
    is.na(hr),
    "",
    paste0(sprintf("%.2f", hr), " (", sprintf("%.2f", low), "-", sprintf("%.2f", high), ")")
  )
}

cohort_day_median <- function(data, variable, day_var = "icu_day") {
  data %>%
    group_by(.data$analysis_cohort, .data[[day_var]]) %>%
    summarise(value = median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(.data$value), NA_real_, .data$value))
}

cohort_median <- function(data, variable) {
  data %>%
    group_by(.data$analysis_cohort) %>%
    summarise(value = median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(.data$value), NA_real_, .data$value))
}

is_variable_usable <- function(data, variable) {
  x <- data[[variable]]
  if (all(is.na(x))) {
    return(FALSE)
  }
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    return(length(unique(x[!is.na(x)])) >= 2)
  }
  length(unique(x[!is.na(x)])) >= 2
}

is_cox_adjustment_usable <- function(data, variable, event_var) {
  if (!is_variable_usable(data, variable)) {
    return(FALSE)
  }
  x <- data[[variable]]
  observed <- x[!is.na(x)]

  if (is.numeric(observed) || is.integer(observed)) {
    unique_values <- sort(unique(observed))
    if (all(unique_values %in% c(0, 1))) {
      counts <- table(observed)
      if (min(counts) < 20) {
        return(FALSE)
      }
      event_balance <- data %>%
        filter(!is.na(.data[[variable]])) %>%
        group_by(.data[[variable]]) %>%
        summarise(
          events = sum(.data[[event_var]], na.rm = TRUE),
          non_events = sum(!.data[[event_var]], na.rm = TRUE),
          .groups = "drop"
        )
      if (any(event_balance$events == 0) || any(event_balance$non_events == 0)) {
        return(FALSE)
      }
    }
  }

  TRUE
}

collapse_rare_factor_levels <- function(x, min_count = 20L) {
  x_chr <- as.character(x)
  tab <- table(x_chr, useNA = "no")
  if (length(tab) == 0) {
    return(factor(x_chr))
  }
  rare_levels <- names(tab)[tab < min_count]
  common_values <- x_chr[!(x_chr %in% rare_levels) & !is.na(x_chr)]
  if (length(common_values) == 0) {
    return(factor(x_chr))
  }
  replacement <- mode_value(common_values)
  x_chr[x_chr %in% rare_levels] <- replacement
  factor(x_chr)
}

# 分析 -----------------------------------------------------------------------

impute_landmark_dynamic_numeric <- function(data, variable, out_variable) {
  missing_variable <- paste0(variable, "_missing_lm")
  within_variable <- paste0(variable, "__within_stay")
  day_median_variable <- paste0(variable, "__cohort_day_median")
  cohort_median_variable <- paste0(variable, "__cohort_median")
  overall_median <- median(data[[variable]], na.rm = TRUE)
  if (is.nan(overall_median)) {
    overall_median <- NA_real_
  }

  day_medians <- cohort_day_median(data, variable, "icu_day") %>%
    rename("{day_median_variable}" := "value")
  cohort_medians <- cohort_median(data, variable) %>%
    rename("{cohort_median_variable}" := "value")

  data %>%
    mutate(
      "{missing_variable}" := as.integer(is.na(.data[[variable]])),
      "{within_variable}" := .data[[variable]]
    ) %>%
    arrange(.data$global_stay_id, .data$icu_day) %>%
    group_by(.data$global_stay_id) %>%
    fill(all_of(within_variable), .direction = "downup") %>%
    ungroup() %>%
    left_join(day_medians, by = c("analysis_cohort", "icu_day")) %>%
    left_join(cohort_medians, by = "analysis_cohort") %>%
    mutate(
      "{out_variable}" := coalesce(
        .data[[within_variable]],
        .data[[day_median_variable]],
        .data[[cohort_median_variable]],
        overall_median
      )
    ) %>%
    select(-all_of(c(within_variable, day_median_variable, cohort_median_variable)))
}

prepare_landmark_subgroup_data <- function() {
  baseline <- read_parquet(file.path(POOLED_DIR, "all_baseline_outcome.parquet")) %>%
    as.data.frame() %>%
    filter(.data$included_final) %>%
    select(
      "global_stay_id",
      baseline_age = "age",
      baseline_sex = "sex",
      baseline_bmi = "bmi",
      baseline_sofa_value = "baseline_sofa",
      baseline_vasopressor_value = "baseline_vasopressor",
      baseline_mv_value = "baseline_mv",
      baseline_copd_value = "copd",
      baseline_ph_value = "baseline_ph",
      baseline_pf_ratio_value = "baseline_pf_ratio",
      baseline_death_28d = "death_28d"
    )

  day_long <- read_parquet(file.path(POOLED_DIR, "all_day_long.parquet")) %>%
    as.data.frame() %>%
    filter(
      .data$included_final,
      .data$landmark_eligible,
      .data$icu_day >= 1L,
      .data$icu_day <= 7L,
      !is.na(.data$twa_paco2)
    ) %>%
    left_join(baseline, by = "global_stay_id") %>%
    mutate(
      landmark_time = pmax(.data$time_from_window_end_to_event_or_censor_days, 1e-06),
      landmark_event = .data$death_28d & .data$time_to_event_28d_days > .data$icu_day,
      analysis_set = unname(cohort_display[.data$analysis_cohort]),
      age_per10 = .data$baseline_age / 10,
      sex_clean = case_when(
        tolower(.data$baseline_sex) %in% c("male", "m") ~ "male",
        tolower(.data$baseline_sex) %in% c("female", "f") ~ "female",
        TRUE ~ "unknown"
      ),
      sex_clean = factor(.data$sex_clean),
      bmi_missing = as.integer(is.na(.data$baseline_bmi)),
      daily_mv_num = as.integer(.data$daily_mv %in% TRUE),
      daily_vasopressor_num = as.integer(.data$daily_vasopressor %in% TRUE),
      daily_crrt_rrt_num = as.integer(.data$daily_crrt_rrt %in% TRUE),
      daily_oxygenation_type_clean = if_else(
        is.na(.data$daily_oxygenation_type) | .data$daily_oxygenation_type == "",
        "unknown",
        .data$daily_oxygenation_type
      ),
      paco2_main_category = case_when(
        .data$twa_paco2 < 35 ~ "<35",
        .data$twa_paco2 <= 50 ~ "35-50",
        TRUE ~ ">50"
      ),
      paco2_main_category = factor(.data$paco2_main_category, levels = c("35-50", "<35", ">50")),
      high_risk_day = as.integer(.data$twa_paco2 < 35 | .data$twa_paco2 > 50),
      subgroup_sex = case_when(
        tolower(.data$baseline_sex) %in% c("male", "m") ~ "Male",
        tolower(.data$baseline_sex) %in% c("female", "f") ~ "Female",
        TRUE ~ NA_character_
      ),
      subgroup_age = if_else(.data$baseline_age < 65, "<65 years", ">=65 years"),
      subgroup_sofa = if_else(.data$baseline_sofa_value < 12, "Low (<12)", "High (>=12)"),
      subgroup_vasopressor = if_else(.data$baseline_vasopressor_value %in% TRUE, "Yes", "No"),
      subgroup_mv = if_else(.data$baseline_mv_value %in% TRUE, "Yes", "No"),
      subgroup_copd = if_else(.data$baseline_copd_value %in% TRUE, "Yes", "No"),
      subgroup_ph = case_when(
        is.na(.data$baseline_ph_value) ~ NA_character_,
        .data$baseline_ph_value < 7.35 ~ "Low (<7.35)",
        TRUE ~ "Non-low (>=7.35)"
      ),
      subgroup_pf = case_when(
        is.na(.data$baseline_pf_ratio_value) ~ NA_character_,
        .data$baseline_pf_ratio_value < 200 ~ "Moderate/severe (<200)",
        TRUE ~ "No or mild (>=200)"
      )
    ) %>%
    arrange(.data$global_stay_id, .data$icu_day) %>%
    group_by(.data$global_stay_id) %>%
    mutate(
      cumulative_high_risk_days = cumsum(.data$high_risk_day),
      any_high_risk = factor(if_else(.data$cumulative_high_risk_days > 0, "Yes", "No"), levels = c("No", "Yes"))
    ) %>%
    ungroup()

  bmi_cohort <- cohort_median(day_long, "baseline_bmi") %>%
    rename(bmi_cohort_median = "value")
  bmi_overall <- median(day_long$baseline_bmi, na.rm = TRUE)
  if (is.nan(bmi_overall)) {
    bmi_overall <- NA_real_
  }

  day_long <- day_long %>%
    left_join(bmi_cohort, by = "analysis_cohort") %>%
    mutate(
      bmi_imp = coalesce(.data$baseline_bmi, .data$bmi_cohort_median, bmi_overall),
      bmi_per5_imp = .data$bmi_imp / 5
    ) %>%
    select(-all_of("bmi_cohort_median"))

  day_long <- impute_landmark_dynamic_numeric(day_long, "daily_sofa", "daily_sofa_imp")
  day_long <- impute_landmark_dynamic_numeric(day_long, "daily_lactate", "daily_lactate_imp")
  day_long <- impute_landmark_dynamic_numeric(day_long, "daily_oxygenation_value", "daily_oxygenation_value_imp")

  day_long %>%
    mutate(
      daily_lactate_imp = pmax(.data$daily_lactate_imp, 0),
      daily_lactate_log1p_imp = log1p(.data$daily_lactate_imp),
      daily_oxygenation_per100_imp = .data$daily_oxygenation_value_imp / 100
    )
}

landmark_adjustment_terms <- function(data, event_var = "landmark_event", exclude_terms = character(0)) {
  model_terms <- c(
    "age_per10",
    "sex_clean",
    "bmi_per5_imp",
    "bmi_missing",
    "daily_sofa_imp",
    "daily_sofa_missing_lm",
    "daily_mv_num",
    "daily_vasopressor_num",
    "daily_crrt_rrt_num",
    "daily_lactate_log1p_imp",
    "daily_lactate_missing_lm",
    "daily_oxygenation_per100_imp",
    "daily_oxygenation_value_missing_lm",
    "daily_oxygenation_type_clean"
  )
  model_terms <- setdiff(model_terms, exclude_terms)
  usable <- model_terms[vapply(model_terms, function(x) is_cox_adjustment_usable(data, x, event_var), logical(1))]
  c(usable, "strata(analysis_cohort)", "strata(icu_day)")
}

patient_counts <- function(data) {
  data %>%
    distinct(.data$global_stay_id, .data$baseline_death_28d) %>%
    summarise(
      n_patients = n(),
      events = sum(.data$baseline_death_28d %in% TRUE),
      .groups = "drop"
    )
}

extract_cox_term <- function(fit, term_pattern) {
  sm <- summary(fit)
  coef_table <- as.data.frame(sm$coefficients)
  ci_table <- as.data.frame(sm$conf.int)
  coef_table$term <- rownames(coef_table)
  ci_table$term <- rownames(ci_table)
  joined <- coef_table %>%
    left_join(ci_table, by = "term") %>%
    filter(grepl(term_pattern, .data$term))
  if (nrow(joined) == 0) {
    return(tibble(adjusted_hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  }
  joined %>%
    slice(1) %>%
    transmute(
      adjusted_hr = exp(.data$coef),
      ci_low = .data$`lower .95`,
      ci_high = .data$`upper .95`,
      p_value = .data$`Pr(>|z|)`
    )
}

subgroup_specs <- tibble(
  subgroup = c(
    "Sex",
    "Age",
    "Baseline SOFA",
    "Baseline vasopressor",
    "Baseline mechanical ventilation",
    "COPD",
    "Baseline pH",
    "PF/oxygenation impairment"
  ),
  variable = c(
    "subgroup_sex",
    "subgroup_age",
    "subgroup_sofa",
    "subgroup_vasopressor",
    "subgroup_mv",
    "subgroup_copd",
    "subgroup_ph",
    "subgroup_pf"
  ),
  levels = I(list(
    c("Male", "Female"),
    c("<65 years", ">=65 years"),
    c("Low (<12)", "High (>=12)"),
    c("No", "Yes"),
    c("No", "Yes"),
    c("No", "Yes"),
    c("Low (<7.35)", "Non-low (>=7.35)"),
    c("No or mild (>=200)", "Moderate/severe (<200)")
  )),
  exclude_terms = I(list(
    "sex_clean",
    "age_per10",
    character(0),
    character(0),
    character(0),
    character(0),
    character(0),
    character(0)
  ))
)

fit_any_high_risk_in_level <- function(data, spec, level) {
  fit_data <- data %>%
    filter(.data[[spec$variable]] == level) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels),
      any_high_risk = factor(.data$any_high_risk, levels = c("No", "Yes"))
    )
  counts <- patient_counts(fit_data)
  if (n_distinct(fit_data$any_high_risk) < 2 || sum(fit_data$landmark_event, na.rm = TRUE) < 10) {
    return(bind_cols(counts, tibble(adjusted_hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)))
  }
  terms <- landmark_adjustment_terms(fit_data, exclude_terms = spec$exclude_terms[[1]])
  form <- as.formula(paste0(
    "Surv(landmark_time, landmark_event) ~ any_high_risk + ",
    paste(c(terms, "cluster(global_stay_id)"), collapse = " + ")
  ))
  fit <- tryCatch(
    suppressWarnings(coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(bind_cols(counts, tibble(adjusted_hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)))
  }
  bind_cols(counts, extract_cox_term(fit, "^any_high_riskYes$"))
}

interaction_p_value <- function(data, spec) {
  fit_data <- data %>%
    filter(!is.na(.data[[spec$variable]])) %>%
    mutate(
      subgroup_model = factor(.data[[spec$variable]], levels = spec$levels[[1]]),
      any_high_risk = factor(.data$any_high_risk, levels = c("No", "Yes")),
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels)
    ) %>%
    filter(!is.na(.data$subgroup_model))
  if (n_distinct(fit_data$subgroup_model) < 2 || n_distinct(fit_data$any_high_risk) < 2) {
    return(NA_real_)
  }
  terms <- landmark_adjustment_terms(fit_data, exclude_terms = spec$exclude_terms[[1]])
  form <- as.formula(paste0(
    "Surv(landmark_time, landmark_event) ~ any_high_risk * subgroup_model + ",
    paste(c(terms, "cluster(global_stay_id)"), collapse = " + ")
  ))
  fit <- tryCatch(
    suppressWarnings(coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NA_real_)
  }
  coef_table <- as.data.frame(summary(fit)$coefficients)
  coef_table$term <- rownames(coef_table)
  interaction_rows <- coef_table %>% filter(grepl("any_high_riskYes:subgroup_model", .data$term, fixed = TRUE))
  if (nrow(interaction_rows) == 0) {
    return(NA_real_)
  }
  as.numeric(interaction_rows$`Pr(>|z|)`[1])
}

fit_subgroup_any_high_risk <- function(data, specs) {
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, ]
    p_int <- interaction_p_value(data, spec)
    bind_rows(lapply(spec$levels[[1]], function(level) {
      est <- fit_any_high_risk_in_level(data, spec, level)
      tibble(
        analysis = "Any cumulative high-risk exposure",
        subgroup = spec$subgroup,
        level = level,
        exposure_contrast = "Any cumulative high-risk exposure vs none",
        n_patients = est$n_patients,
        events = est$events,
        adjusted_hr = est$adjusted_hr,
        ci_low = est$ci_low,
        ci_high = est$ci_high,
        p_value = est$p_value,
        p_interaction = p_int
      )
    }))
  }))
}

impute_pamm_dynamic_numeric <- function(data, variable, out_variable) {
  missing_variable <- paste0(variable, "_missing")
  within_variable <- paste0(variable, "__within_stay")
  day_median_variable <- paste0(variable, "__cohort_day_median")
  cohort_median_variable <- paste0(variable, "__cohort_median")
  overall_median <- median(data[[variable]], na.rm = TRUE)
  if (is.nan(overall_median)) {
    overall_median <- NA_real_
  }

  day_medians <- data %>%
    group_by(.data$analysis_cohort, .data$source_icu_day) %>%
    summarise("{day_median_variable}" := median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate("{day_median_variable}" := if_else(is.nan(.data[[day_median_variable]]), NA_real_, .data[[day_median_variable]]))

  cohort_medians <- cohort_median(data, variable) %>%
    rename("{cohort_median_variable}" := "value")

  data %>%
    mutate(
      "{missing_variable}" := as.integer(is.na(.data[[variable]])),
      "{within_variable}" := .data[[variable]]
    ) %>%
    arrange(.data$global_stay_id, .data$ped_interval) %>%
    group_by(.data$global_stay_id) %>%
    fill(all_of(within_variable), .direction = "downup") %>%
    ungroup() %>%
    left_join(day_medians, by = c("analysis_cohort", "source_icu_day")) %>%
    left_join(cohort_medians, by = "analysis_cohort") %>%
    mutate(
      "{out_variable}" := coalesce(
        .data[[within_variable]],
        .data[[day_median_variable]],
        .data[[cohort_median_variable]],
        overall_median
      )
    ) %>%
    select(-all_of(c(within_variable, day_median_variable, cohort_median_variable)))
}

prepare_pamm_data <- function(data, scope_id) {
  fit_data <- data
  if (scope_id != "pooled") {
    fit_data <- fit_data %>% filter(.data$analysis_cohort == scope_to_cohort[[scope_id]])
  }

  bmi_cohort <- cohort_median(fit_data, "bmi") %>%
    rename(bmi_cohort_median = "value")
  bmi_overall <- median(fit_data$bmi, na.rm = TRUE)
  if (is.nan(bmi_overall)) {
    bmi_overall <- NA_real_
  }

  fit_data <- fit_data %>%
    left_join(bmi_cohort, by = "analysis_cohort") %>%
    mutate(
      scope_id = scope_id,
      sex_clean = case_when(
        tolower(.data$sex) %in% c("male", "m") ~ "male",
        tolower(.data$sex) %in% c("female", "f") ~ "female",
        TRUE ~ "unknown"
      ),
      sex_clean = factor(.data$sex_clean),
      bmi_missing = as.integer(is.na(.data$bmi)),
      bmi_imp = coalesce(.data$bmi, .data$bmi_cohort_median, bmi_overall),
      age_per10 = .data$age / 10,
      bmi_per5_imp = .data$bmi_imp / 5,
      daily_mv_num = as.integer(.data$daily_mv %in% TRUE),
      daily_vasopressor_num = as.integer(.data$daily_vasopressor %in% TRUE),
      daily_crrt_rrt_num = as.integer(.data$daily_crrt_rrt %in% TRUE),
      daily_oxygenation_type_clean = if_else(
        is.na(.data$daily_oxygenation_type) | .data$daily_oxygenation_type == "",
        "unknown",
        .data$daily_oxygenation_type
      ),
      daily_oxygenation_type_clean = factor(.data$daily_oxygenation_type_clean),
      baseline_stratum = if_else(
        .data$scope_id == "pooled",
        paste(.data$analysis_cohort, .data$ped_interval, sep = "__day"),
        paste0("day", .data$ped_interval)
      ),
      baseline_stratum = factor(.data$baseline_stratum)
    ) %>%
    select(-"bmi_cohort_median")

  fit_data <- impute_pamm_dynamic_numeric(fit_data, "daily_sofa", "daily_sofa_imp")
  fit_data <- impute_pamm_dynamic_numeric(fit_data, "daily_lactate", "daily_lactate_imp")
  fit_data <- impute_pamm_dynamic_numeric(fit_data, "daily_oxygenation_value", "daily_oxygenation_value_imp")

  ph_center <- median(fit_data$daily_ph, na.rm = TRUE)
  fit_data %>%
    mutate(
      daily_lactate_imp = pmax(.data$daily_lactate_imp, 0),
      daily_lactate_log1p_imp = log1p(.data$daily_lactate_imp),
      daily_oxygenation_per100_imp = .data$daily_oxygenation_value_imp / 100,
      daily_ph_centered = .data$daily_ph - ph_center
    )
}

pamm_model_terms <- function(ph_adjusted = FALSE) {
  terms <- c(
    "baseline_stratum",
    "age_per10",
    "sex_clean",
    "bmi_per5_imp",
    "bmi_missing",
    "daily_sofa_imp",
    "daily_sofa_missing",
    "daily_mv_num",
    "daily_vasopressor_num",
    "daily_crrt_rrt_num",
    "daily_lactate_log1p_imp",
    "daily_lactate_missing",
    "daily_oxygenation_per100_imp",
    "daily_oxygenation_value_missing",
    "daily_oxygenation_type_clean"
  )
  if (ph_adjusted) {
    terms <- c(terms, "daily_ph_centered")
  }
  terms
}

pamm_formula <- function(data, ph_adjusted = FALSE, paco2_k = PRIMARY_K) {
  terms <- pamm_model_terms(ph_adjusted)
  usable_terms <- terms[vapply(terms, function(x) is_variable_usable(data, x), logical(1))]
  exposure_term <- paste0("s(twa_paco2, k = ", paco2_k, ", bs = 'cr')")
  as.formula(paste0(
    "ped_status ~ ",
    paste(c(usable_terms, exposure_term), collapse = " + "),
    " + offset(offset_log_interval)"
  ))
}

pamm_template_newdata <- function(fit_data, grid) {
  template <- tibble(
    twa_paco2 = grid,
    source_icu_day = 4,
    offset_log_interval = 0,
    age_per10 = median(fit_data$age_per10, na.rm = TRUE),
    bmi_per5_imp = median(fit_data$bmi_per5_imp, na.rm = TRUE),
    bmi_missing = as.integer(round(mean(fit_data$bmi_missing, na.rm = TRUE))),
    daily_sofa_imp = median(fit_data$daily_sofa_imp, na.rm = TRUE),
    daily_sofa_missing = as.integer(round(mean(fit_data$daily_sofa_missing, na.rm = TRUE))),
    daily_mv_num = as.integer(as.numeric(mode_value(fit_data$daily_mv_num))),
    daily_vasopressor_num = as.integer(as.numeric(mode_value(fit_data$daily_vasopressor_num))),
    daily_crrt_rrt_num = as.integer(as.numeric(mode_value(fit_data$daily_crrt_rrt_num))),
    daily_lactate_log1p_imp = median(fit_data$daily_lactate_log1p_imp, na.rm = TRUE),
    daily_lactate_missing = as.integer(round(mean(fit_data$daily_lactate_missing, na.rm = TRUE))),
    daily_oxygenation_per100_imp = median(fit_data$daily_oxygenation_per100_imp, na.rm = TRUE),
    daily_oxygenation_value_missing = as.integer(round(mean(fit_data$daily_oxygenation_value_missing, na.rm = TRUE))),
    daily_ph_centered = 0
  )

  template$sex_clean <- factor(mode_value(fit_data$sex_clean), levels = levels(fit_data$sex_clean))
  template$daily_oxygenation_type_clean <- factor(
    mode_value(fit_data$daily_oxygenation_type_clean),
    levels = levels(fit_data$daily_oxygenation_type_clean)
  )
  template$baseline_stratum <- factor(
    mode_value(fit_data$baseline_stratum),
    levels = levels(fit_data$baseline_stratum)
  )
  template
}

predict_pamm_curve <- function(fit, fit_data, model_label, scope_id, ph_adjusted, grid_n = 281) {
  p1 <- as.numeric(quantile(fit_data$twa_paco2, 0.01, na.rm = TRUE))
  p99 <- as.numeric(quantile(fit_data$twa_paco2, 0.99, na.rm = TRUE))
  grid_min <- max(20, floor(p1))
  grid_max <- min(90, ceiling(p99))
  if (grid_max <= grid_min) {
    grid_min <- floor(min(fit_data$twa_paco2, na.rm = TRUE))
    grid_max <- ceiling(max(fit_data$twa_paco2, na.rm = TRUE))
  }

  grid <- seq(grid_min, grid_max, length.out = grid_n)
  nd <- pamm_template_newdata(fit_data, grid)
  xmat <- predict(fit, nd, type = "lpmatrix")
  eta <- as.vector(xmat %*% coef(fit))
  nadir_idx <- which.min(eta)
  xref <- xmat[nadir_idx, , drop = FALSE]
  xdiff <- sweep(xmat, 2, xref[1, ], "-")
  vc <- vcov(fit)
  se <- sqrt(pmax(0, rowSums((xdiff %*% vc) * xdiff)))
  log_hr <- eta - eta[nadir_idx]

  tibble(
    model = model_label,
    scope_id = scope_id,
    analysis_set = unname(scope_display[scope_id]),
    ph_adjusted = ph_adjusted,
    paco2 = grid,
    support_p1 = p1,
    support_p99 = p99,
    nadir_paco2 = grid[nadir_idx],
    log_hr_vs_nadir = log_hr,
    se_log_hr_vs_nadir = se,
    hr_vs_nadir = exp(log_hr),
    hr_low_vs_nadir = exp(log_hr - 1.96 * se),
    hr_high_vs_nadir = exp(log_hr + 1.96 * se)
  )
}

nearest_curve_point <- function(curve, point) {
  curve %>%
    slice_min(abs(.data$paco2 - point), n = 1, with_ties = FALSE) %>%
    transmute(
      "{paste0('hr_at_', point, '_vs_nadir')}" := .data$hr_vs_nadir,
      "{paste0('hr_low_at_', point, '_vs_nadir')}" := .data$hr_low_vs_nadir,
      "{paste0('hr_high_at_', point, '_vs_nadir')}" := .data$hr_high_vs_nadir
    )
}

summarise_pamm_curve <- function(curve, fit_data, fit, form) {
  near5 <- curve %>% filter(.data$hr_vs_nadir <= 1.05)
  bind_cols(
    tibble(
      model = first(curve$model),
      analysis_set = first(curve$analysis_set),
      n_patients = n_distinct(fit_data$global_stay_id),
      n_ped_rows = nrow(fit_data),
      events = sum(fit_data$ped_status, na.rm = TRUE),
      nadir_paco2 = first(curve$nadir_paco2),
      near_min_5pct_low = if (nrow(near5) > 0) min(near5$paco2) else NA_real_,
      near_min_5pct_high = if (nrow(near5) > 0) max(near5$paco2) else NA_real_,
      aic = AIC(fit),
      model_formula = paste(deparse(form), collapse = " ")
    ),
    nearest_curve_point(curve, 35),
    nearest_curve_point(curve, 50),
    nearest_curve_point(curve, 60)
  )
}

fit_pamm_sensitivity <- function(base_data, scope_id, ph_adjusted) {
  fit_data <- prepare_pamm_data(base_data, scope_id)
  model_label <- if (ph_adjusted) "pH-adjusted Model 3" else "Primary pH-available Model 3"
  form <- pamm_formula(fit_data, ph_adjusted = ph_adjusted, paco2_k = PRIMARY_K)
  message(
    "Fitting ", model_label, " / ", scope_display[[scope_id]],
    " (rows=", nrow(fit_data),
    ", stays=", n_distinct(fit_data$global_stay_id),
    ", events=", sum(fit_data$ped_status, na.rm = TRUE), ")"
  )

  fit <- bam(
    form,
    data = fit_data,
    family = poisson(link = "log"),
    method = "fREML",
    discrete = TRUE,
    nthreads = ANALYSIS_NTHREADS
  )
  curve <- predict_pamm_curve(fit, fit_data, model_label, scope_id, ph_adjusted)
  summary <- summarise_pamm_curve(curve, fit_data, fit, form)
  saveRDS(
    fit,
    file.path(MODEL_OUT_DIR, paste0("result4_", ifelse(ph_adjusted, "ph_adjusted", "primary_ph_available"), "_", scope_id, ".rds"))
  )
  list(curve = curve, summary = summary)
}

landmark_data <- prepare_landmark_subgroup_data()
subgroup_main_raw <- fit_subgroup_any_high_risk(landmark_data, subgroup_specs)
subgroup_all_raw <- subgroup_main_raw

ph_available_ped <- read_parquet(file.path(PAMM_DIR, "ped_sensitivity_ph.parquet")) %>%
  as.data.frame() %>%
  filter(!is.na(.data$twa_paco2), !is.na(.data$daily_ph))

pamm_fit_results <- list()
fit_index <- 1L
for (ph_adjusted in c(FALSE, TRUE)) {
  pamm_fit_results[[fit_index]] <- fit_pamm_sensitivity(ph_available_ped, "pooled", ph_adjusted)
  fit_index <- fit_index + 1L
}

ph_curves <- bind_rows(lapply(pamm_fit_results, `[[`, "curve")) %>%
  mutate(
    model = factor(.data$model, levels = c("Primary pH-available Model 3", "pH-adjusted Model 3")),
    analysis_set = factor(.data$analysis_set, levels = unname(scope_display))
  )

ph_summaries <- bind_rows(lapply(pamm_fit_results, `[[`, "summary")) %>%
  mutate(
    model = factor(.data$model, levels = c("Primary pH-available Model 3", "pH-adjusted Model 3")),
    analysis_set = factor(.data$analysis_set, levels = unname(scope_display))
  ) %>%
  arrange(.data$analysis_set, .data$model)

# 导出结果 -------------------------------------------------------------------

table_s6 <- subgroup_all_raw %>%
  mutate(
    Subgroup = .data$subgroup,
    Level = .data$level,
    `N patients` = .data$n_patients,
    Events = .data$events,
    `Adjusted HR` = if_else(is.na(.data$adjusted_hr), "", sprintf("%.2f", .data$adjusted_hr)),
    `95% CI` = fmt_hr_ci_english(.data$adjusted_hr, .data$ci_low, .data$ci_high) %>%
      gsub("^.*\\((.*)\\)$", "\\1", .),
    `P value` = fmt_p(.data$p_value),
    `P for interaction` = if_else(!duplicated(.data$Subgroup), fmt_p(.data$p_interaction), "")
  ) %>%
  select(
    "Subgroup",
    "Level",
    "N patients",
    "Events",
    "Adjusted HR",
    "95% CI",
    "P value",
    "P for interaction"
  )

table_s7 <- ph_summaries %>%
  mutate(
    Model = as.character(.data$model),
    `N patients` = .data$n_patients,
    `N PED rows` = .data$n_ped_rows,
    Events = .data$events,
    `Nadir PaCO2, mmHg` = sprintf("%.1f", .data$nadir_paco2),
    `5% near-minimum range, mmHg` = paste0(sprintf("%.1f", .data$near_min_5pct_low), "-", sprintf("%.1f", .data$near_min_5pct_high)),
    `HR at 35 mmHg` = fmt_hr_ci_english(.data$hr_at_35_vs_nadir, .data$hr_low_at_35_vs_nadir, .data$hr_high_at_35_vs_nadir),
    `HR at 50 mmHg` = fmt_hr_ci_english(.data$hr_at_50_vs_nadir, .data$hr_low_at_50_vs_nadir, .data$hr_high_at_50_vs_nadir),
    `HR at 60 mmHg` = fmt_hr_ci_english(.data$hr_at_60_vs_nadir, .data$hr_low_at_60_vs_nadir, .data$hr_high_at_60_vs_nadir),
    AIC = round(.data$aic, 1)
  ) %>%
  select(
    "Model",
    "N patients",
    "N PED rows",
    "Events",
    "Nadir PaCO2, mmHg",
    "5% near-minimum range, mmHg",
    "HR at 35 mmHg",
    "HR at 50 mmHg",
    "HR at 60 mmHg",
    "AIC"
  )

figure3_data <- subgroup_main_raw %>%
  mutate(
    subgroup = factor(.data$subgroup, levels = subgroup_specs$subgroup),
    row_label = paste0(.data$subgroup, ": ", .data$level),
    row_label = factor(.data$row_label, levels = rev(paste0(
      rep(subgroup_specs$subgroup, lengths(subgroup_specs$levels)),
      ": ",
      unlist(subgroup_specs$levels)
    ))),
    hr_label = if_else(
      is.na(.data$adjusted_hr),
      "",
      paste0(sprintf("%.2f", .data$adjusted_hr), " (", sprintf("%.2f", .data$ci_low), "-", sprintf("%.2f", .data$ci_high), ")")
    ),
    p_int_label = if_else(
      !duplicated(.data$subgroup),
      if_else(is.na(.data$p_interaction), "", paste0("P int = ", fmt_p_plain(.data$p_interaction))),
      ""
    )
  )

figure3 <- ggplot(figure3_data, aes(x = .data$adjusted_hr, y = .data$row_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_errorbar(aes(xmin = .data$ci_low, xmax = .data$ci_high), orientation = "y", width = 0, linewidth = 0.55, color = "#2C7FB8", na.rm = TRUE) +
  geom_point(size = 2.0, color = "#2C7FB8", na.rm = TRUE) +
  geom_text(aes(x = 3.25, label = .data$hr_label), hjust = 0, size = 2.75, color = "#222222", na.rm = TRUE) +
  geom_text(aes(x = 9.2, label = .data$p_int_label), hjust = 0, size = 2.75, color = "#222222", na.rm = TRUE) +
  scale_x_log10(limits = c(0.55, 13), breaks = c(0.75, 1, 1.5, 2, 3)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Adjusted hazard ratio for any cumulative high-risk PaCO2 exposure",
    y = NULL
  ) +
  base_theme(9) +
  theme(
    plot.margin = margin(5.5, 120, 5.5, 45),
    legend.position = "none"
  )

range_background <- function() {
  list(
    annotate("rect", xmin = 35, xmax = 50, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.55),
    annotate("rect", xmin = 40, xmax = 45, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.45),
    geom_vline(xintercept = c(35, 50), linetype = "dashed", color = "grey40", linewidth = 0.35),
    geom_vline(xintercept = c(40, 45), linetype = "dotted", color = "grey40", linewidth = 0.35)
  )
}

figure_s6_data <- ph_curves %>%
  filter(.data$analysis_set == "Pooled")

figure_s6 <- ggplot(figure_s6_data, aes(x = .data$paco2, y = .data$hr_vs_nadir, color = .data$model, fill = .data$model)) +
  range_background() +
  geom_ribbon(aes(ymin = .data$hr_low_vs_nadir, ymax = .data$hr_high_vs_nadir), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  scale_color_manual(values = figure_colors) +
  scale_fill_manual(values = figure_colors) +
  coord_cartesian(xlim = c(20, 80), ylim = c(0.75, 3.2)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "Adjusted hazard ratio vs model-specific nadir"
  ) +
  base_theme(10)

ggsave(
  file.path(FIG_MAIN_DIR, "Figure_3_subgroup_analysis.png"),
  figure3,
  width = 10.5,
  height = 6.8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_MAIN_DIR, "Figure_3_subgroup_analysis.pdf"),
  figure3,
  width = 10.5,
  height = 6.8,
  bg = "white"
)

ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S6_ph_adjusted_sensitivity_pamm.png"),
  figure_s6,
  width = 7.5,
  height = 4.8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S6_ph_adjusted_sensitivity_pamm.pdf"),
  figure_s6,
  width = 7.5,
  height = 4.8,
  bg = "white"
)

write_csv(table_s6, file.path(TABLE_SUPP_DIR, "Table_S7_subgroup_estimates.csv"), na = "")
write.xlsx(table_s6, file.path(TABLE_SUPP_DIR, "Table_S7_subgroup_estimates.xlsx"), overwrite = TRUE)
write_csv(subgroup_all_raw, file.path(TABLE_SUPP_DIR, "Result_4_subgroup_estimates_raw.csv"), na = "")

write_csv(table_s7, file.path(TABLE_SUPP_DIR, "Table_S8_ph_adjusted_sensitivity.csv"), na = "")
write.xlsx(table_s7, file.path(TABLE_SUPP_DIR, "Table_S8_ph_adjusted_sensitivity.xlsx"), overwrite = TRUE)
write_csv(ph_curves, file.path(TABLE_SUPP_DIR, "Figure_S6_ph_adjusted_sensitivity_pamm_curves.csv"), na = "")
write_csv(ph_summaries, file.path(TABLE_SUPP_DIR, "Table_S8_ph_adjusted_sensitivity_raw.csv"), na = "")

main_subgroups <- subgroup_main_raw %>%
  filter(!is.na(.data$adjusted_hr))

interaction_ps <- subgroup_main_raw %>%
  distinct(.data$subgroup, .data$p_interaction)

copd_yes <- subgroup_main_raw %>% filter(.data$subgroup == "COPD", .data$level == "Yes") %>% slice(1)
mv_yes <- subgroup_main_raw %>% filter(.data$subgroup == "Baseline mechanical ventilation", .data$level == "Yes") %>% slice(1)
ph_low <- subgroup_main_raw %>% filter(.data$subgroup == "Baseline pH", .data$level == "Low (<7.35)") %>% slice(1)
pf_impair <- subgroup_main_raw %>% filter(.data$subgroup == "PF/oxygenation impairment", .data$level == "Moderate/severe (<200)") %>% slice(1)

pooled_primary_ph <- ph_summaries %>%
  filter(.data$analysis_set == "Pooled", .data$model == "Primary pH-available Model 3") %>%
  slice(1)
pooled_ph_adjusted <- ph_summaries %>%
  filter(.data$analysis_set == "Pooled", .data$model == "pH-adjusted Model 3") %>%
  slice(1)

result4_text <- c(
  "## Result 4. 亚组分析与 pH 调整敏感性分析",
  "",
  paste0(
    "预设亚组分析按性别、年龄、基线 SOFA、基线血管活性药物使用、基线机械通气、COPD、基线 pH 和 PF/氧合损伤分层。总体上，任意累计高风险 PaCO2 暴露与 28 天全因死亡的关联在各亚组中的点估计均位于 HR>1，但部分置信区间跨越 1；各亚组 adjusted HR 点估计范围为 ",
    sprintf("%.2f", min(main_subgroups$adjusted_hr, na.rm = TRUE)),
    "-",
    sprintf("%.2f", max(main_subgroups$adjusted_hr, na.rm = TRUE)),
    "。各预设亚组均未观察到统计学明确的交互证据（Figure 3；Table S7）。在 COPD 患者中，任意累计高风险 PaCO2 暴露的 adjusted HR 为 ",
    fmt_estimate(copd_yes$adjusted_hr, copd_yes$ci_low, copd_yes$ci_high),
    "；在基线机械通气患者中为 ",
    fmt_estimate(mv_yes$adjusted_hr, mv_yes$ci_low, mv_yes$ci_high),
    "；在基线 pH <7.35 患者中为 ",
    fmt_estimate(ph_low$adjusted_hr, ph_low$ci_low, ph_low$ci_high),
    "；在 PF <200 患者中为 ",
    fmt_estimate(pf_impair$adjusted_hr, pf_impair$ci_low, pf_impair$ci_high),
    "（Figure 3；Table S7）。"
  ),
  "",
  paste0(
    "pH-adjusted 敏感性 PAMM 中纳入 pH 可用的 ",
    fmt_int(pooled_ph_adjusted$n_patients),
    " 例患者，对应 ",
    fmt_int(pooled_ph_adjusted$n_ped_rows),
    " 个 PED 行和 ",
    fmt_int(pooled_ph_adjusted$events),
    " 个事件。额外加入 pH 后最低风险锚点为 ",
    sprintf("%.1f", pooled_ph_adjusted$nadir_paco2),
    " mmHg，5% 近最低风险区间为 ",
    sprintf("%.1f", pooled_ph_adjusted$near_min_5pct_low),
    "-",
    sprintf("%.1f", pooled_ph_adjusted$near_min_5pct_high),
    " mmHg。pH-adjusted Model 3 中，TWA-PaCO2 为 35、50 和 60 mmHg 时相对于模型最低点的 adjusted HR 分别为 ",
    fmt_estimate(pooled_ph_adjusted$hr_at_35_vs_nadir, pooled_ph_adjusted$hr_low_at_35_vs_nadir, pooled_ph_adjusted$hr_high_at_35_vs_nadir),
    "、",
    fmt_estimate(pooled_ph_adjusted$hr_at_50_vs_nadir, pooled_ph_adjusted$hr_low_at_50_vs_nadir, pooled_ph_adjusted$hr_high_at_50_vs_nadir),
    "和 ",
    fmt_estimate(pooled_ph_adjusted$hr_at_60_vs_nadir, pooled_ph_adjusted$hr_low_at_60_vs_nadir, pooled_ph_adjusted$hr_high_at_60_vs_nadir),
    "。加入 pH 后高 PaCO2 端关联减弱，但未观察到较高 PaCO2 与更低死亡风险相关（Figure S6；Table S8）。"
  )
)

writeLines(
  result4_text,
  con = file.path(MANUSCRIPT_DIR, "simulated_results_result4_zh.md"),
  useBytes = TRUE
)

result4_qc <- tibble(
  item = c(
    "subgroup_landmark_rows",
    "subgroup_landmark_patients",
    "ph_available_ped_rows",
    "ph_available_patients",
    "pooled_ph_adjusted_nadir_paco2",
    "pooled_ph_adjusted_hr_at_50_vs_nadir"
  ),
  value = c(
    nrow(landmark_data),
    n_distinct(landmark_data$global_stay_id),
    nrow(ph_available_ped),
    n_distinct(ph_available_ped$global_stay_id),
    round(pooled_ph_adjusted$nadir_paco2, 2),
    round(pooled_ph_adjusted$hr_at_50_vs_nadir, 3)
  )
)
write_csv(result4_qc, file.path(TABLE_SUPP_DIR, "Result_4_qc_summary.csv"), na = "")

message("Done Result 4.")
message("Figure 3: ", file.path(FIG_MAIN_DIR, "Figure_3_subgroup_analysis.png"))
message("Figure S6: ", file.path(FIG_SUPP_DIR, "Figure_S6_ph_adjusted_sensitivity_pamm.png"))
message("Table S7: ", file.path(TABLE_SUPP_DIR, "Table_S7_subgroup_estimates.xlsx"))
message("Table S8: ", file.path(TABLE_SUPP_DIR, "Table_S8_ph_adjusted_sensitivity.xlsx"))
message("Result text: ", file.path(MANUSCRIPT_DIR, "simulated_results_result4_zh.md"))
