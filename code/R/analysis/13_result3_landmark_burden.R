options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(readr)
  library(survival)
  library(tibble)
  library(tidyr)
})

# Configuration --------------------------------------------------------------

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

FIG_MAIN_DIR <- file.path(RESULT_ROOT, "figures", "main")
FIG_SUPP_DIR <- file.path(RESULT_ROOT, "figures", "supplementary")
TABLE_MAIN_DIR <- file.path(RESULT_ROOT, "tables", "main")
TABLE_SUPP_DIR <- file.path(RESULT_ROOT, "tables", "supplementary")
MANUSCRIPT_DIR <- file.path(RESULT_ROOT, "manuscript")

dir.create(FIG_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

NOMINAL_LANDMARK_WINDOW_MINUTES <- 24 * 60

message("Running Result 3: landmark categories and high-risk PaCO2 burden...")

cohort_display <- c(
  MIMIC = "MIMIC-IV",
  AmsterdamUMCdb = "AmsterdamUMCdb",
  `Chinese cohort` = "Chinese"
)

analysis_set_order <- c("Pooled", "MIMIC-IV", "AmsterdamUMCdb", "Chinese")
scope_ids <- c("pooled", "mimic", "amsterdam", "chinese")
scope_display <- c(
  pooled = "Pooled",
  mimic = "MIMIC-IV",
  amsterdam = "AmsterdamUMCdb",
  chinese = "Chinese"
)
scope_to_cohort <- c(
  mimic = "MIMIC",
  amsterdam = "AmsterdamUMCdb",
  chinese = "Chinese cohort"
)

landmark_model_display <- c(
  model1 = "Model 1",
  model2 = "Model 2",
  model3 = "Model 3"
)

figure_colors <- c(
  "Pooled" = "#222222",
  "MIMIC-IV" = "#2C7FB8",
  "AmsterdamUMCdb" = "#41AB5D",
  "Chinese" = "#D95F0E"
)

base_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.tag = element_text(face = "bold", size = base_size + 3)
    )
}

# Shared helpers -------------------------------------------------------------

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  names(sort(table(x), decreasing = TRUE))[1]
}

cohort_day_median <- function(data, variable) {
  data %>%
    group_by(.data$analysis_cohort, .data$icu_day) %>%
    summarise(value = median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(.data$value), NA_real_, .data$value))
}

cohort_median <- function(data, variable) {
  data %>%
    group_by(.data$analysis_cohort) %>%
    summarise(value = median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(.data$value), NA_real_, .data$value))
}

impute_dynamic_numeric <- function(data, variable, out_variable) {
  missing_variable <- paste0(variable, "_missing_lm")
  within_variable <- paste0(variable, "__within_stay")
  day_median_variable <- paste0(variable, "__cohort_day_median")
  cohort_median_variable <- paste0(variable, "__cohort_median")
  overall_median <- median(data[[variable]], na.rm = TRUE)
  if (is.nan(overall_median)) {
    overall_median <- NA_real_
  }

  day_medians <- cohort_day_median(data, variable) %>%
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
    fill(all_of(within_variable), .direction = "down") %>%
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

impute_baseline_numeric <- function(data, variable, out_variable, missing_variable = paste0(variable, "_missing")) {
  cohort_medians <- cohort_median(data, variable) %>%
    rename("{variable}__cohort_median" := "value")
  cohort_var <- paste0(variable, "__cohort_median")
  overall_median <- median(data[[variable]], na.rm = TRUE)
  if (is.nan(overall_median)) {
    overall_median <- NA_real_
  }

  data %>%
    left_join(cohort_medians, by = "analysis_cohort") %>%
    mutate(
      "{missing_variable}" := as.integer(is.na(.data[[variable]])),
      "{out_variable}" := coalesce(.data[[variable]], .data[[cohort_var]], overall_median)
    ) %>%
    select(-all_of(cohort_var))
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

filter_scope <- function(data, scope_id) {
  if (scope_id == "pooled") {
    return(data)
  }
  data %>% filter(.data$analysis_cohort == scope_to_cohort[[scope_id]])
}

format_p <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "P < .001",
    TRUE ~ paste0("P = ", sub("^0", "", sprintf("%.3f", x)))
  )
}

format_num <- function(x, digits = 2) {
  if_else(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_estimate <- function(hr, low, high) {
  paste0(sprintf("%.2f", hr), "（95% CI ", sprintf("%.2f", low), "-", sprintf("%.2f", high), "）")
}

# Landmark category data -----------------------------------------------------

prepare_landmark_data <- function() {
  baseline <- read_parquet(file.path(POOLED_DIR, "all_baseline_outcome.parquet")) %>%
    as.data.frame() %>%
    filter(.data$included_final) %>%
    select(
      "global_stay_id",
      baseline_age = "age",
      baseline_sex = "sex",
      baseline_bmi = "bmi"
    )

  day_long_source <- read_parquet(file.path(POOLED_DIR, "all_day_long.parquet")) %>%
    as.data.frame() %>%
    filter(
      .data$included_final,
      .data$icu_day >= 1L,
      .data$icu_day <= 7L,
      !is.na(.data$twa_paco2)
    )

  landmark_risk_set_qc <- day_long_source %>%
    mutate(
      complete_nominal_window = near(.data$window_duration_min, NOMINAL_LANDMARK_WINDOW_MINUTES),
      strict_landmark_eligible = .data$landmark_eligible & .data$complete_nominal_window
    ) %>%
    group_by(.data$analysis_cohort, .data$icu_day) %>%
    summarise(
      previously_eligible_rows = sum(.data$landmark_eligible, na.rm = TRUE),
      strict_landmark_rows = sum(.data$strict_landmark_eligible, na.rm = TRUE),
      excluded_incomplete_window_rows = sum(.data$landmark_eligible & !.data$complete_nominal_window, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(
    landmark_risk_set_qc,
    file.path(TABLE_SUPP_DIR, "Result_3_fixed_landmark_risk_set_qc.csv"),
    na = ""
  )

  day_long <- day_long_source %>%
    filter(
      .data$landmark_eligible,
      near(.data$window_duration_min, NOMINAL_LANDMARK_WINDOW_MINUTES)
    ) %>%
    left_join(baseline, by = "global_stay_id") %>%
    mutate(
      landmark_time = pmax(.data$time_from_window_end_to_event_or_censor_days, 1e-06),
      landmark_event = .data$death_28d & .data$time_from_window_end_to_event_or_censor_days > 0,
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
      paco2_fine_category = case_when(
        .data$twa_paco2 < 35 ~ "<35",
        .data$twa_paco2 <= 40 ~ "35-40",
        .data$twa_paco2 <= 45 ~ "40-45",
        .data$twa_paco2 <= 50 ~ "45-50",
        TRUE ~ ">50"
      ),
      paco2_fine_category = factor(
        .data$paco2_fine_category,
        levels = c("40-45", "<35", "35-40", "45-50", ">50")
      )
    )

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

  day_long <- impute_dynamic_numeric(day_long, "daily_sofa", "daily_sofa_imp")
  day_long <- impute_dynamic_numeric(day_long, "daily_lactate", "daily_lactate_imp")
  day_long <- impute_dynamic_numeric(day_long, "daily_oxygenation_value", "daily_oxygenation_value_imp")

  day_long %>%
    mutate(
      daily_lactate_imp = pmax(.data$daily_lactate_imp, 0),
      daily_lactate_log1p_imp = log1p(.data$daily_lactate_imp),
      daily_oxygenation_per100_imp = .data$daily_oxygenation_value_imp / 100
    )
}

landmark_adjustment_terms <- function(data, scope_id, model_id = "model3", include_day_strata = TRUE) {
  model_terms <- switch(
    model_id,
    model1 = character(0),
    model2 = c(
      "age_per10",
      "sex_clean",
      "bmi_per5_imp",
      "bmi_missing"
    ),
    model3 = c(
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
    ),
    stop("Unsupported landmark model_id: ", model_id, call. = FALSE)
  )

  usable <- model_terms[vapply(model_terms, function(x) is_cox_adjustment_usable(data, x, "landmark_event"), logical(1))]
  if (scope_id == "pooled") {
    usable <- c(usable, "strata(analysis_cohort)")
  }
  if (include_day_strata) {
    usable <- c(usable, "strata(icu_day)")
  }
  usable
}

landmark_covariate_label <- function(model_id) {
  switch(
    model_id,
    model1 = "Exposure + cohort/day strata",
    model2 = "Model 1 + age + sex + BMI",
    model3 = paste(
      "Model 2 + daily SOFA + daily mechanical ventilation + daily vasopressor use",
      "+ daily CRRT/RRT + daily lactate + daily oxygenation status"
    ),
    stop("Unsupported landmark model_id: ", model_id, call. = FALSE)
  )
}

fit_landmark_main <- function(data, scope_id, model_id = "model3") {
  model_label <- unname(landmark_model_display[[model_id]])
  covariate_label <- landmark_covariate_label(model_id)

  fit_data <- data %>%
    filter_scope(scope_id) %>%
    filter(!is.na(.data$paco2_main_category)) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels),
      exposure_group = droplevels(.data$paco2_main_category)
    )

  rhs <- paste(
    c("exposure_group", landmark_adjustment_terms(fit_data, scope_id, model_id), "cluster(global_stay_id)"),
    collapse = " + "
  )
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))
  fit <- coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50))
  sm <- summary(fit)
  coef_table <- as.data.frame(sm$coefficients)
  ci_table <- as.data.frame(sm$conf.int)
  coef_table$term <- rownames(coef_table)
  ci_table$term <- rownames(ci_table)
  se_col <- intersect(c("robust se", "se(coef)"), colnames(coef_table))[1]

  coef_table %>%
    left_join(ci_table, by = "term") %>%
    filter(grepl("^exposure_group", .data$term)) %>%
    transmute(
      scope_id = .env$scope_id,
      model_id = .env$model_id,
      model = .env$model_label,
      analysis_set = unname(scope_display[scope_id]),
      time_window = "Days 1-7",
      n_days_estimated = n_distinct(fit_data$icu_day),
      exposure_category = sub("^exposure_group", "", .data$term),
      reference = "35-50",
      n_rows = nrow(fit_data),
      n_patients = n_distinct(fit_data$global_stay_id),
      events = sum(fit_data$landmark_event, na.rm = TRUE),
      log_hr = .data$coef,
      se_log_hr = .data[[se_col]],
      adjusted_hr = exp(.data$coef),
      ci_low = .data$`lower .95`,
      ci_high = .data$`upper .95`,
      p_value = .data$`Pr(>|z|)`,
      aic = AIC(fit),
      covariates = .env$covariate_label,
      model_formula = paste(deparse(form), collapse = " ")
    )
}

extract_cox_term_rows <- function(fit, fit_data, term_map, id_columns, count_variable = NULL) {
  expected <- tibble(
    term = names(term_map),
    category_or_contrast = unname(term_map)
  )

  if (is.null(fit)) {
    estimates <- expected %>%
      transmute(
        category_or_contrast,
        log_hr = NA_real_,
        se_log_hr = NA_real_,
        adjusted_hr = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        aic = NA_real_
      )
  } else {
    sm <- summary(fit)
    coef_table <- as.data.frame(sm$coefficients)
    ci_table <- as.data.frame(sm$conf.int)
    coef_table$term <- rownames(coef_table)
    ci_table$term <- rownames(ci_table)
    se_col <- intersect(c("robust se", "se(coef)"), colnames(coef_table))[1]

    estimates <- coef_table %>%
      left_join(ci_table, by = "term") %>%
      filter(.data$term %in% names(term_map)) %>%
      transmute(
        category_or_contrast = unname(term_map[.data$term]),
        log_hr = .data$coef,
        se_log_hr = .data[[se_col]],
        adjusted_hr = exp(.data$coef),
        ci_low = .data$`lower .95`,
        ci_high = .data$`upper .95`,
        p_value = .data$`Pr(>|z|)`,
        aic = AIC(fit)
      )

    estimates <- expected %>%
      select("category_or_contrast") %>%
      left_join(estimates, by = "category_or_contrast")
  }

  if (!is.null(count_variable) && count_variable %in% names(fit_data)) {
    counts <- fit_data %>%
      mutate(category_or_contrast = as.character(.data[[count_variable]])) %>%
      group_by(.data$category_or_contrast) %>%
      summarise(
        n_in_category = n_distinct(.data$global_stay_id),
        events_in_category = sum(.data$landmark_event, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    counts <- tibble(
      category_or_contrast = unname(term_map),
      n_in_category = n_distinct(fit_data$global_stay_id),
      events_in_category = sum(fit_data$landmark_event, na.rm = TRUE)
    )
  }

  bind_cols(id_columns[rep(1, nrow(expected)), , drop = FALSE], expected %>% select("category_or_contrast")) %>%
    left_join(estimates, by = "category_or_contrast") %>%
    left_join(counts, by = "category_or_contrast") %>%
    mutate(
      n_at_landmark = n_distinct(fit_data$global_stay_id),
      events_at_landmark = sum(fit_data$landmark_event, na.rm = TRUE)
    )
}

fit_landmark_day <- function(data, scope_id, day, model_id = "model3") {
  model_label <- unname(landmark_model_display[[model_id]])
  covariate_label <- landmark_covariate_label(model_id)

  fit_data <- data %>%
    filter_scope(scope_id) %>%
    filter(.data$icu_day == .env$day, !is.na(.data$paco2_main_category)) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels),
      exposure_group = droplevels(.data$paco2_main_category)
    )

  rhs <- paste(
    c("exposure_group", landmark_adjustment_terms(fit_data, scope_id, model_id, include_day_strata = FALSE)),
    collapse = " + "
  )
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))

  fit <- if (nlevels(fit_data$exposure_group) >= 2 && length(unique(fit_data$landmark_event)) >= 2) {
    tryCatch(
      coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50)),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  extract_cox_term_rows(
    fit,
    fit_data,
    c(exposure_group35 = "<35", exposure_group50 = ">50") %>%
      setNames(c("exposure_group<35", "exposure_group>50")),
    tibble(
      scope_id = scope_id,
      model_id = model_id,
      model = model_label,
      analysis_set = unname(scope_display[scope_id]),
      icu_day = day,
      time_window = paste0("Day ", day),
      reference = "35-50",
      covariates = covariate_label,
      model_formula = if (is.null(fit)) paste(deparse(form), collapse = " ") else paste(deparse(form), collapse = " ")
    ),
    count_variable = "exposure_group"
  ) %>%
    rename(exposure_category = "category_or_contrast") %>%
    mutate(
      n_rows = .data$n_at_landmark,
      n_patients = .data$n_at_landmark,
      events = .data$events_at_landmark
    )
}

fit_fine_landmark_day <- function(data, scope_id, day, model_id = "model3") {
  model_label <- unname(landmark_model_display[[model_id]])
  covariate_label <- landmark_covariate_label(model_id)

  fit_data <- data %>%
    filter_scope(scope_id) %>%
    filter(.data$icu_day == .env$day, !is.na(.data$paco2_fine_category)) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels),
      exposure_group = droplevels(.data$paco2_fine_category)
    )

  rhs <- paste(
    c("exposure_group", landmark_adjustment_terms(fit_data, scope_id, model_id, include_day_strata = FALSE)),
    collapse = " + "
  )
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))

  fit <- if (nlevels(fit_data$exposure_group) >= 2 && length(unique(fit_data$landmark_event)) >= 2) {
    tryCatch(
      coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50)),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  extract_cox_term_rows(
    fit,
    fit_data,
    c(
      `exposure_group<35` = "<35",
      `exposure_group35-40` = "35-40",
      `exposure_group45-50` = "45-50",
      `exposure_group>50` = ">50"
    ),
    tibble(
      scope_id = scope_id,
      model_id = model_id,
      model = model_label,
      analysis_set = unname(scope_display[scope_id]),
      icu_day = day,
      time_window = paste0("Day ", day),
      reference = "40-45",
      covariates = covariate_label,
      model_formula = paste(deparse(form), collapse = " ")
    ),
    count_variable = "exposure_group"
  ) %>%
    rename(exposure_category = "category_or_contrast") %>%
    mutate(
      n_rows = .data$n_at_landmark,
      n_patients = .data$n_at_landmark,
      events = .data$events_at_landmark
    )
}

fit_fine_landmark_stacked <- function(data, scope_id, model_id = "model3") {
  model_label <- unname(landmark_model_display[[model_id]])
  covariate_label <- landmark_covariate_label(model_id)

  fit_data <- data %>%
    filter_scope(scope_id) %>%
    filter(!is.na(.data$paco2_fine_category)) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels),
      exposure_group = droplevels(.data$paco2_fine_category)
    )

  rhs <- paste(
    c("exposure_group", landmark_adjustment_terms(fit_data, scope_id, model_id), "cluster(global_stay_id)"),
    collapse = " + "
  )
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))
  fit <- coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50))

  extract_cox_term_rows(
    fit,
    fit_data,
    c(
      `exposure_group<35` = "<35",
      `exposure_group35-40` = "35-40",
      `exposure_group45-50` = "45-50",
      `exposure_group>50` = ">50"
    ),
    tibble(
      scope_id = scope_id,
      model_id = model_id,
      model = model_label,
      analysis_set = unname(scope_display[scope_id]),
      icu_day = NA_integer_,
      time_window = "Days 1-7",
      n_days_estimated = n_distinct(fit_data$icu_day),
      reference = "40-45",
      covariates = covariate_label,
      model_formula = paste(deparse(form), collapse = " ")
    ),
    count_variable = "exposure_group"
  ) %>%
    rename(exposure_category = "category_or_contrast") %>%
    mutate(
      n_rows = nrow(fit_data),
      n_patients = n_distinct(fit_data$global_stay_id),
      events = sum(fit_data$landmark_event, na.rm = TRUE)
    )
}

# Patient-level burden data --------------------------------------------------

prepare_burden_data <- function() {
  baseline <- read_parquet(file.path(POOLED_DIR, "all_baseline_outcome.parquet")) %>%
    as.data.frame() %>%
    filter(.data$included_final) %>%
    mutate(
      analysis_set = unname(cohort_display[.data$analysis_cohort]),
      outcome_time = pmax(.data$time_to_event_28d_days, 1e-06),
      outcome_event = .data$death_28d %in% TRUE,
      age_per10 = .data$age / 10,
      sex_clean = case_when(
        tolower(.data$sex) %in% c("male", "m") ~ "male",
        tolower(.data$sex) %in% c("female", "f") ~ "female",
        TRUE ~ "unknown"
      ),
      sex_clean = factor(.data$sex_clean),
      baseline_mv_num = as.integer(.data$baseline_mv %in% TRUE),
      baseline_vasopressor_num = as.integer(.data$baseline_vasopressor %in% TRUE),
      baseline_crrt_rrt_num = as.integer(.data$baseline_crrt_rrt %in% TRUE),
      baseline_oxygenation_type_clean = case_when(
        !is.na(.data$baseline_pf_ratio) ~ "pf_ratio",
        !is.na(.data$baseline_pao2) ~ "pao2",
        TRUE ~ "unknown"
      ),
      baseline_oxygenation_type_clean = factor(.data$baseline_oxygenation_type_clean),
      baseline_oxygenation_value = coalesce(.data$baseline_pf_ratio, .data$baseline_pao2)
    )

  baseline <- impute_baseline_numeric(baseline, "bmi", "bmi_imp", "bmi_missing") %>%
    mutate(bmi_per5_imp = .data$bmi_imp / 5)
  baseline <- impute_baseline_numeric(baseline, "baseline_sofa", "baseline_sofa_imp", "baseline_sofa_missing")
  baseline <- impute_baseline_numeric(baseline, "baseline_lactate", "baseline_lactate_imp", "baseline_lactate_missing")
  baseline <- impute_baseline_numeric(
    baseline,
    "baseline_oxygenation_value",
    "baseline_oxygenation_value_imp",
    "baseline_oxygenation_value_missing"
  )

  burden <- read_parquet(file.path(POOLED_DIR, "all_day_long.parquet")) %>%
    as.data.frame() %>%
    filter(
      .data$included_final,
      .data$observable_icu_day,
      .data$icu_day >= 1L,
      .data$icu_day <= 7L,
      !is.na(.data$twa_paco2)
    ) %>%
    mutate(
      high_risk_day = .data$twa_paco2 < 35 | .data$twa_paco2 > 50,
      deviation_dose = pmax(35 - .data$twa_paco2, 0) + pmax(.data$twa_paco2 - 50, 0)
    ) %>%
    group_by(.data$global_stay_id) %>%
    summarise(
      n_observable_days = n(),
      high_risk_days = sum(.data$high_risk_day, na.rm = TRUE),
      cumulative_deviation_dose = sum(.data$deviation_dose, na.rm = TRUE),
      .groups = "drop"
    )

  baseline %>%
    left_join(burden, by = "global_stay_id") %>%
    mutate(
      n_observable_days = coalesce(.data$n_observable_days, 0L),
      high_risk_days = coalesce(.data$high_risk_days, 0L),
      cumulative_deviation_dose = coalesce(.data$cumulative_deviation_dose, 0),
      any_high_risk = if_else(.data$high_risk_days > 0, "Yes", "No"),
      any_high_risk = factor(.data$any_high_risk, levels = c("No", "Yes")),
      high_risk_days_cat = case_when(
        .data$high_risk_days == 0 ~ "0",
        .data$high_risk_days == 1 ~ "1",
        .data$high_risk_days <= 3 ~ "2-3",
        TRUE ~ ">=4"
      ),
      high_risk_days_cat = factor(.data$high_risk_days_cat, levels = c("0", "1", "2-3", ">=4")),
      high_risk_prop = if_else(.data$n_observable_days > 0, .data$high_risk_days / .data$n_observable_days, NA_real_),
      high_risk_prop_cat = case_when(
        .data$high_risk_prop == 0 ~ "0%",
        .data$high_risk_prop <= 0.25 ~ ">0-25%",
        .data$high_risk_prop <= 0.50 ~ ">25-50%",
        TRUE ~ ">50%"
      ),
      high_risk_prop_cat = factor(.data$high_risk_prop_cat, levels = c("0%", ">0-25%", ">25-50%", ">50%")),
      cumulative_deviation_per10 = .data$cumulative_deviation_dose / 10,
      baseline_lactate_imp = pmax(.data$baseline_lactate_imp, 0),
      baseline_lactate_log1p_imp = log1p(.data$baseline_lactate_imp),
      baseline_oxygenation_per100_imp = .data$baseline_oxygenation_value_imp / 100
    )
}

burden_adjustment_terms <- function(data, scope_id) {
  terms <- c(
    "age_per10",
    "sex_clean",
    "bmi_per5_imp",
    "bmi_missing",
    "baseline_sofa_imp",
    "baseline_sofa_missing",
    "baseline_mv_num",
    "baseline_vasopressor_num",
    "baseline_crrt_rrt_num",
    "baseline_lactate_log1p_imp",
    "baseline_lactate_missing",
    "baseline_oxygenation_per100_imp",
    "baseline_oxygenation_value_missing",
    "baseline_oxygenation_type_clean"
  )
  usable <- terms[vapply(terms, function(x) is_cox_adjustment_usable(data, x, "outcome_event"), logical(1))]
  if (scope_id == "pooled") {
    usable <- c(usable, "strata(analysis_cohort)")
  }
  usable
}

fit_burden_model <- function(data, exposure_term, scope_id = "pooled") {
  fit_data <- data %>%
    filter_scope(scope_id) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels)
    )

  rhs <- paste(c(exposure_term, landmark_adjustment_terms(fit_data, scope_id), "cluster(global_stay_id)"), collapse = " + ")
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))
  fit <- coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50))
  list(fit = fit, data = fit_data, formula = form)
}

extract_burden_rows <- function(fit_obj, exposure_metric, category_map, reference, trend_p = NA_real_) {
  sm <- summary(fit_obj$fit)
  coef_table <- as.data.frame(sm$coefficients)
  ci_table <- as.data.frame(sm$conf.int)
  coef_table$term <- rownames(coef_table)
  ci_table$term <- rownames(ci_table)

  estimate_rows <- coef_table %>%
    left_join(ci_table, by = "term") %>%
    filter(.data$term %in% names(category_map)) %>%
    transmute(
      exposure_metric = exposure_metric,
      category_or_contrast = unname(category_map[.data$term]),
      reference = reference,
      n_rows = nrow(fit_obj$data),
      n_patients = n_distinct(fit_obj$data$global_stay_id),
      events = sum(fit_obj$data$landmark_event, na.rm = TRUE),
      adjusted_hr = exp(.data$coef),
      ci_low = .data$`lower .95`,
      ci_high = .data$`upper .95`,
      p_value = .data$`Pr(>|z|)`,
      p_for_trend = trend_p,
      model_formula = paste(deparse(fit_obj$formula), collapse = " ")
    )

  estimate_rows
}

trend_p_value <- function(data, exposure_term, scope_id = "pooled") {
  fit_obj <- fit_burden_model(data, exposure_term, scope_id)
  coef_table <- as.data.frame(summary(fit_obj$fit)$coefficients)
  coef_table[1, "Pr(>|z|)"]
}

category_summary <- function(data, variable, levels, metric_label) {
  data %>%
    mutate(category = factor(.data[[variable]], levels = levels)) %>%
    group_by(.data$category) %>%
    summarise(
      exposure_metric = metric_label,
      category_or_contrast = as.character(first(.data$category)),
      n_rows_category = n(),
      n_patients_category = n_distinct(.data$global_stay_id),
      events_category = sum(.data$landmark_event, na.rm = TRUE),
      .groups = "drop"
    )
}

category_summary_by_day <- function(data, variable, levels, metric_label) {
  day_totals <- data %>%
    group_by(.data$icu_day) %>%
    summarise(n_at_landmark = n_distinct(.data$global_stay_id), .groups = "drop")

  data %>%
    mutate(category_or_contrast = factor(.data[[variable]], levels = levels)) %>%
    group_by(.data$icu_day, .data$category_or_contrast) %>%
    summarise(
      exposure_metric = metric_label,
      n_in_category = n_distinct(.data$global_stay_id),
      events_in_category = sum(.data$landmark_event, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(day_totals, by = "icu_day") %>%
    mutate(category_or_contrast = as.character(.data$category_or_contrast))
}

category_distribution_by_day <- function(data, variable, levels, scheme_label) {
  count_data <- bind_rows(
    data %>% mutate(analysis_set_for_count = "Pooled"),
    data %>% mutate(analysis_set_for_count = .data$analysis_set)
  )
  day_totals <- count_data %>%
    group_by(.data$analysis_set_for_count, .data$icu_day) %>%
    summarise(n_at_landmark = n_distinct(.data$global_stay_id), .groups = "drop")

  count_data %>%
    mutate(exposure_category = factor(.data[[variable]], levels = levels)) %>%
    group_by(.data$analysis_set_for_count, .data$icu_day, .data$exposure_category) %>%
    summarise(
      n_patients = n_distinct(.data$global_stay_id),
      events = sum(.data$landmark_event, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(day_totals, by = c("analysis_set_for_count", "icu_day")) %>%
    mutate(
      Scheme = scheme_label,
      `Analysis set` = factor(.data$analysis_set_for_count, levels = analysis_set_order),
      `ICU day` = paste0("Day ", .data$icu_day),
      `Exposure category` = as.character(.data$exposure_category),
      `N at landmark` = .data$n_at_landmark,
      `N in category` = .data$n_patients,
      Events = .data$events,
      `Event rate` = if_else(.data$n_patients > 0, paste0(sprintf("%.1f", .data$events / .data$n_patients * 100), "%"), "")
    ) %>%
    arrange(.data$`Analysis set`, .data$icu_day, .data$exposure_category) %>%
    transmute(
      Scheme,
      `Analysis set` = as.character(.data$`Analysis set`),
      `ICU day`,
      `Exposure category`,
      `N at landmark`,
      `N in category`,
      Events,
      `Event rate`
    )
}

fit_burden_day_model <- function(data, exposure_term, scope_id = "pooled", day) {
  fit_data <- data %>%
    filter_scope(scope_id) %>%
    filter(.data$icu_day == .env$day) %>%
    mutate(
      daily_oxygenation_type_clean = collapse_rare_factor_levels(.data$daily_oxygenation_type_clean),
      across(where(is.factor), droplevels)
    )

  rhs <- paste(
    c(exposure_term, landmark_adjustment_terms(fit_data, scope_id, "model3", include_day_strata = FALSE)),
    collapse = " + "
  )
  form <- as.formula(paste0("Surv(landmark_time, landmark_event) ~ ", rhs))

  fit <- if (is_variable_usable(fit_data, exposure_term) && length(unique(fit_data$landmark_event)) >= 2) {
    tryCatch(
      coxph(form, data = fit_data, ties = "efron", control = coxph.control(iter.max = 50)),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  list(fit = fit, data = fit_data, formula = form)
}

extract_burden_day_rows <- function(data, exposure_term, exposure_metric, category_map, reference, scope_id = "pooled", day) {
  fit_obj <- fit_burden_day_model(data, exposure_term, scope_id, day)

  extract_cox_term_rows(
    fit_obj$fit,
    fit_obj$data,
    category_map,
    tibble(
      scope_id = scope_id,
      analysis_set = unname(scope_display[scope_id]),
      icu_day = day,
      time_window = paste0("Day ", day),
      exposure_metric = exposure_metric,
      reference = reference,
      model_formula = paste(deparse(fit_obj$formula), collapse = " ")
    ),
    count_variable = if (is.factor(fit_obj$data[[exposure_term]]) || is.character(fit_obj$data[[exposure_term]])) {
      exposure_term
    } else {
      NULL
    }
  )
}

range_label <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return("")
  }
  if (length(unique(x)) == 1) {
    return(fmt_int(unique(x)))
  }
  paste0(fmt_int(min(x)), "-", fmt_int(max(x)))
}

estimate_label <- function(hr, low, high) {
  if (is.na(hr) || is.na(low) || is.na(high)) {
    return("")
  }
  paste0(format_num(hr, 2), " (", format_num(low, 2), "-", format_num(high, 2), ")")
}

# Run analyses ---------------------------------------------------------------

landmark_data <- prepare_landmark_data()
stacked_landmark_hierarchy_raw <- bind_rows(lapply(names(landmark_model_display), function(model_id) {
  fit_landmark_main(landmark_data, "pooled", model_id)
})) %>%
  mutate(
    model = factor(.data$model, levels = unname(landmark_model_display)),
    exposure_category = factor(.data$exposure_category, levels = c("<35", ">50"))
  ) %>%
  arrange(.data$model, .data$exposure_category) %>%
  mutate(
    model = as.character(.data$model),
    exposure_category = as.character(.data$exposure_category)
  )

stacked_cohort_hierarchy_raw <- bind_rows(lapply(setdiff(scope_ids, "pooled"), function(scope_id) {
  bind_rows(lapply(names(landmark_model_display), function(model_id) {
    fit_landmark_main(landmark_data, scope_id, model_id)
  }))
})) %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order),
    model = factor(.data$model, levels = unname(landmark_model_display)),
    exposure_category = factor(.data$exposure_category, levels = c("<35", ">50"))
  ) %>%
  arrange(.data$analysis_set, .data$model, .data$exposure_category) %>%
  mutate(
    analysis_set = as.character(.data$analysis_set),
    model = as.character(.data$model),
    exposure_category = as.character(.data$exposure_category)
  )

landmark_day_raw <- bind_rows(lapply(scope_ids, function(scope_id) {
  bind_rows(lapply(1:7, function(day) {
    fit_landmark_day(landmark_data, scope_id, day, "model3")
  }))
})) %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order),
    exposure_category = factor(.data$exposure_category, levels = c("<35", ">50"))
  ) %>%
  arrange(.data$analysis_set, .data$exposure_category) %>%
  mutate(
    analysis_set = as.character(.data$analysis_set),
    exposure_category = as.character(.data$exposure_category)
  )

landmark_raw <- bind_rows(stacked_landmark_hierarchy_raw, stacked_cohort_hierarchy_raw) %>%
  filter(.data$model_id == "model3") %>%
  mutate(
    icu_day = NA_integer_,
    summary_method = "Stacked landmark Cox model with patient-clustered robust standard errors"
  ) %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order),
    exposure_category = factor(.data$exposure_category, levels = c("<35", ">50"))
  ) %>%
  arrange(.data$analysis_set, .data$exposure_category) %>%
  mutate(
    analysis_set = as.character(.data$analysis_set),
    exposure_category = as.character(.data$exposure_category)
  )

landmark_detail_raw <- bind_rows(
  landmark_day_raw %>% mutate(row_type = "Day-specific"),
  landmark_raw %>% mutate(row_type = "Stacked overall")
)

fine_day_raw <- bind_rows(lapply(scope_ids, function(scope_id) {
  bind_rows(lapply(1:7, function(day) {
    fit_fine_landmark_day(landmark_data, scope_id, day, "model3")
  }))
})) %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order),
    exposure_category = factor(.data$exposure_category, levels = c("<35", "35-40", "45-50", ">50"))
  ) %>%
  arrange(.data$analysis_set, .data$exposure_category) %>%
  mutate(
    analysis_set = as.character(.data$analysis_set),
    exposure_category = as.character(.data$exposure_category)
  )

fine_raw <- bind_rows(lapply(scope_ids, function(scope_id) {
  fit_fine_landmark_stacked(landmark_data, scope_id, "model3")
})) %>%
  mutate(summary_method = "Stacked landmark Cox model with patient-clustered robust standard errors") %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order),
    exposure_category = factor(.data$exposure_category, levels = c("<35", "35-40", "45-50", ">50"))
  ) %>%
  arrange(.data$analysis_set, .data$exposure_category) %>%
  mutate(
    analysis_set = as.character(.data$analysis_set),
    exposure_category = as.character(.data$exposure_category)
  )

fine_detail_raw <- bind_rows(
  fine_day_raw %>% mutate(row_type = "Day-specific"),
  fine_raw %>% mutate(row_type = "Stacked overall")
)

burden_landmark <- landmark_data %>%
  arrange(.data$global_stay_id, .data$icu_day) %>%
  group_by(.data$global_stay_id) %>%
  mutate(
    high_risk_day = .data$twa_paco2 < 35 | .data$twa_paco2 > 50,
    deviation_dose = pmax(35 - .data$twa_paco2, 0) + pmax(.data$twa_paco2 - 50, 0),
    cumulative_high_risk_days = cumsum(.data$high_risk_day),
    cumulative_deviation_dose = cumsum(.data$deviation_dose),
    cumulative_high_risk_prop = .data$cumulative_high_risk_days / row_number(),
    any_high_risk = if_else(.data$cumulative_high_risk_days > 0, "Yes", "No"),
    any_high_risk = factor(.data$any_high_risk, levels = c("No", "Yes")),
    high_risk_days_cat = case_when(
      .data$cumulative_high_risk_days == 0 ~ "0",
      .data$cumulative_high_risk_days == 1 ~ "1",
      .data$cumulative_high_risk_days <= 3 ~ "2-3",
      TRUE ~ ">=4"
    ),
    high_risk_days_cat = factor(.data$high_risk_days_cat, levels = c("0", "1", "2-3", ">=4")),
    high_risk_prop_cat = case_when(
      .data$cumulative_high_risk_prop == 0 ~ "0%",
      .data$cumulative_high_risk_prop <= 0.25 ~ ">0-25%",
      .data$cumulative_high_risk_prop <= 0.50 ~ ">25-50%",
      TRUE ~ ">50%"
    ),
    high_risk_prop_cat = factor(.data$high_risk_prop_cat, levels = c("0%", ">0-25%", ">25-50%", ">50%")),
    cumulative_deviation_per10 = .data$cumulative_deviation_dose / 10
  ) %>%
  ungroup()

burden_pooled <- filter_scope(burden_landmark, "pooled")

ever_high_risk_summary <- burden_pooled %>%
  group_by(.data$global_stay_id) %>%
  summarise(ever_high_risk = any(.data$high_risk_day %in% TRUE), .groups = "drop") %>%
  summarise(
    n_patients = n(),
    any_exposed_patients = sum(.data$ever_high_risk, na.rm = TRUE),
    any_exposed_pct = .data$any_exposed_patients / .data$n_patients * 100,
    .groups = "drop"
  )

burden_day_raw <- bind_rows(lapply(1:7, function(day) {
  bind_rows(
    extract_burden_day_rows(
      burden_landmark,
      "any_high_risk",
      "Any high-risk exposure",
      c(any_high_riskYes = "Yes"),
      "No",
      "pooled",
      day
    ),
    extract_burden_day_rows(
      burden_landmark,
      "high_risk_days_cat",
      "Number of high-risk days",
      c(high_risk_days_cat1 = "1", `high_risk_days_cat2-3` = "2-3", `high_risk_days_cat>=4` = ">=4"),
      "0",
      "pooled",
      day
    ),
    extract_burden_day_rows(
      burden_landmark,
      "high_risk_prop_cat",
      "Proportion of high-risk days",
      c(
        `high_risk_prop_cat>0-25%` = ">0-25%",
        `high_risk_prop_cat>25-50%` = ">25-50%",
        `high_risk_prop_cat>50%` = ">50%"
      ),
      "0%",
      "pooled",
      day
    ),
    extract_burden_day_rows(
      burden_landmark,
      "cumulative_deviation_per10",
      "Cumulative deviation dose",
      c(cumulative_deviation_per10 = "per 10 mmHg-days"),
      "0 mmHg-days",
      "pooled",
      day
    )
  )
}))

days_trend_p <- trend_p_value(burden_landmark, "cumulative_high_risk_days", "pooled")
prop_trend_p <- trend_p_value(burden_landmark, "cumulative_high_risk_prop", "pooled")

burden_raw <- bind_rows(
  extract_burden_rows(
    fit_burden_model(burden_landmark, "any_high_risk", "pooled"),
    "Any high-risk exposure",
    c(any_high_riskYes = "Yes"),
    "No"
  ),
  extract_burden_rows(
    fit_burden_model(burden_landmark, "high_risk_days_cat", "pooled"),
    "Number of high-risk days",
    c(high_risk_days_cat1 = "1", `high_risk_days_cat2-3` = "2-3", `high_risk_days_cat>=4` = ">=4"),
    "0",
    trend_p = days_trend_p
  ),
  extract_burden_rows(
    fit_burden_model(burden_landmark, "high_risk_prop_cat", "pooled"),
    "Proportion of high-risk days",
    c(
      `high_risk_prop_cat>0-25%` = ">0-25%",
      `high_risk_prop_cat>25-50%` = ">25-50%",
      `high_risk_prop_cat>50%` = ">50%"
    ),
    "0%",
    trend_p = prop_trend_p
  ),
  extract_burden_rows(
    fit_burden_model(burden_landmark, "cumulative_deviation_per10", "pooled"),
    "Cumulative deviation dose",
    c(cumulative_deviation_per10 = "per 10 mmHg-days"),
    "0 mmHg-days"
  )
) %>%
  mutate(
    scope_id = "pooled",
    analysis_set = "Pooled",
    icu_day = NA_integer_,
    time_window = "Days 1-7",
    n_days_estimated = n_distinct(burden_landmark$icu_day),
    log_hr = log(.data$adjusted_hr),
    se_log_hr = (log(.data$ci_high) - log(.data$ci_low)) / (2 * qnorm(0.975)),
    summary_method = "Stacked landmark Cox model with patient-clustered robust standard errors"
  )

burden_raw <- burden_raw %>%
  mutate(
    p_for_trend = case_when(
      .data$exposure_metric == "Number of high-risk days" ~ .env$days_trend_p,
      .data$exposure_metric == "Proportion of high-risk days" ~ .env$prop_trend_p,
      TRUE ~ NA_real_
    )
  )

burden_counts_by_day <- bind_rows(
  category_summary_by_day(burden_pooled, "any_high_risk", c("No", "Yes"), "Any high-risk exposure"),
  category_summary_by_day(burden_pooled, "high_risk_days_cat", c("0", "1", "2-3", ">=4"), "Number of high-risk days"),
  category_summary_by_day(
    burden_pooled,
    "high_risk_prop_cat",
    c("0%", ">0-25%", ">25-50%", ">50%"),
    "Proportion of high-risk days"
  ),
  burden_pooled %>%
    group_by(.data$icu_day) %>%
    summarise(
      exposure_metric = "Cumulative deviation dose",
      category_or_contrast = "per 10 mmHg-days",
      n_at_landmark = n_distinct(.data$global_stay_id),
      n_in_category = n_distinct(.data$global_stay_id),
      events_in_category = sum(.data$landmark_event, na.rm = TRUE),
      .groups = "drop"
    )
)

burden_counts_summary <- burden_counts_by_day %>%
  group_by(.data$exposure_metric, .data$category_or_contrast) %>%
  summarise(
    n_at_landmark_range = range_label(.data$n_at_landmark),
    n_in_category_range = range_label(.data$n_in_category),
    events_in_category_range = range_label(.data$events_in_category),
    .groups = "drop"
  )

burden_raw <- burden_raw %>%
  left_join(burden_counts_summary, by = c("exposure_metric", "category_or_contrast"))

# Formatted tables -----------------------------------------------------------

landmark_count_summary <- landmark_day_raw %>%
  group_by(.data$analysis_set, .data$exposure_category, .data$reference) %>%
  summarise(
    n_at_landmark_range = range_label(.data$n_at_landmark),
    n_in_category_range = range_label(.data$n_in_category),
    events_in_category_range = range_label(.data$events_in_category),
    .groups = "drop"
  )

table2 <- landmark_raw %>%
  left_join(landmark_count_summary, by = c("analysis_set", "exposure_category", "reference")) %>%
  transmute(
    `Analysis set` = .data$analysis_set,
    `Exposure category` = .data$exposure_category,
    `N at landmark range` = .data$n_at_landmark_range,
    `N in category range` = .data$n_in_category_range,
    `Events in category range` = .data$events_in_category_range,
    `Adjusted HR` = format_num(.data$adjusted_hr, 2),
    `95% CI` = if_else(is.na(.data$ci_low), "", paste0(format_num(.data$ci_low, 2), "-", format_num(.data$ci_high, 2))),
    `P value` = format_p(.data$p_value)
  )

fine_count_summary <- fine_day_raw %>%
  group_by(.data$analysis_set, .data$exposure_category, .data$reference) %>%
  summarise(
    n_at_landmark_range = range_label(.data$n_at_landmark),
    n_in_category_range = range_label(.data$n_in_category),
    events_in_category_range = range_label(.data$events_in_category),
    .groups = "drop"
  )

table_s3 <- fine_raw %>%
  left_join(fine_count_summary, by = c("analysis_set", "exposure_category", "reference")) %>%
  transmute(
    `Analysis set` = .data$analysis_set,
    `Exposure category` = .data$exposure_category,
    `N at landmark range` = .data$n_at_landmark_range,
    `N in category range` = .data$n_in_category_range,
    `Events in category range` = .data$events_in_category_range,
    `Adjusted HR` = format_num(.data$adjusted_hr, 2),
    `95% CI` = if_else(is.na(.data$ci_low), "", paste0(format_num(.data$ci_low, 2), "-", format_num(.data$ci_high, 2))),
    `P value` = format_p(.data$p_value)
  )

table_s4 <- bind_rows(
  category_distribution_by_day(landmark_data, "paco2_main_category", c("<35", "35-50", ">50"), "Main category"),
  category_distribution_by_day(landmark_data, "paco2_fine_category", c("<35", "35-40", "40-45", "45-50", ">50"), "Fine category")
)

reference_rows <- burden_counts_summary %>%
  filter(
    (.data$exposure_metric == "Any high-risk exposure" & .data$category_or_contrast == "No") |
      (.data$exposure_metric == "Number of high-risk days" & .data$category_or_contrast == "0") |
      (.data$exposure_metric == "Proportion of high-risk days" & .data$category_or_contrast == "0%")
  ) %>%
  transmute(
    scope_id = "pooled",
    analysis_set = "Pooled",
    exposure_metric,
    category_or_contrast,
    reference = "Reference",
    time_window = "Days 1-7",
    n_days_estimated = 7L,
    n_at_landmark_range,
    n_in_category_range,
    events_in_category_range,
    adjusted_hr = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = NA_real_,
    p_for_trend = NA_real_,
    log_hr = NA_real_,
    se_log_hr = NA_real_,
    model_formula = ""
  )

table3_raw <- bind_rows(reference_rows, burden_raw)

table3_template <- read_csv(
  file.path(SCRIPT_ROOT, "assets", "templates", "Table_3_high_risk_paco2_burden_template.csv"),
  show_col_types = FALSE
) %>%
  select("Exposure metric", "Category or contrast", "Reference")

table3 <- table3_template %>%
  left_join(
    table3_raw,
    by = c(
      "Exposure metric" = "exposure_metric",
      "Category or contrast" = "category_or_contrast",
      "Reference" = "reference"
    )
  ) %>%
  transmute(
    `Exposure metric`,
    `Category or contrast`,
    Reference,
    `Landmark days` = "Day 1-7",
    `N at landmark range` = .data$n_at_landmark_range,
    `N in category range` = .data$n_in_category_range,
    `Events in category range` = .data$events_in_category_range,
    `N days estimated` = .data$n_days_estimated,
    `Adjusted HR` = if_else(is.na(.data$adjusted_hr), "Reference", format_num(.data$adjusted_hr, 2)),
    `95% CI` = if_else(is.na(.data$ci_low), "", paste0(format_num(.data$ci_low, 2), "-", format_num(.data$ci_high, 2))),
    `P value` = format_p(.data$p_value),
    `P for trend` = format_p(.data$p_for_trend)
  )

table_s5 <- read_csv(
  file.path(SCRIPT_ROOT, "assets", "templates", "Table_S5_secondary_exposure_definitions_template.csv"),
  show_col_types = FALSE
)

# Figures --------------------------------------------------------------------

exposure_colors <- c(
  "<35" = "#2C7FB8",
  ">50" = "#D95F0E"
)

figure2_data <- bind_rows(
  landmark_day_raw %>%
    filter(.data$analysis_set == "Pooled") %>%
    mutate(day_label = paste0("Day ", .data$icu_day), point_type = "Day-specific"),
  landmark_raw %>%
    filter(.data$analysis_set == "Pooled") %>%
    mutate(day_label = "Overall", point_type = "Stacked overall")
) %>%
  mutate(
    day_label = factor(.data$day_label, levels = rev(c(paste0("Day ", 1:7), "Overall"))),
    exposure_category = factor(.data$exposure_category, levels = c("<35", ">50")),
    hr_label = if_else(
      is.na(.data$adjusted_hr),
      "",
      paste0(sprintf("%.2f", .data$adjusted_hr), " (", sprintf("%.2f", .data$ci_low), "-", sprintf("%.2f", .data$ci_high), ")")
    )
  )

figure2 <- ggplot(figure2_data, aes(x = .data$adjusted_hr, y = .data$day_label, color = .data$exposure_category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_errorbar(aes(xmin = .data$ci_low, xmax = .data$ci_high), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(aes(shape = .data$point_type), size = 2.2) +
  geom_text(aes(x = 2.35, label = .data$hr_label), hjust = 0, size = 2.7, color = "#222222") +
  facet_wrap(~exposure_category, ncol = 2) +
  scale_color_manual(values = exposure_colors) +
  scale_shape_manual(values = c("Day-specific" = 16, "Stacked overall" = 18)) +
  scale_x_log10(limits = c(0.55, 5.5), breaks = c(0.75, 1, 1.5, 2, 3)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Adjusted hazard ratio vs 35-50 mmHg",
    y = NULL
  ) +
  base_theme(10) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5))

figure_s5_data <- landmark_raw %>%
  filter(.data$analysis_set != "Pooled") %>%
  mutate(
    analysis_set = factor(.data$analysis_set, levels = analysis_set_order[-1]),
    label = paste(.data$analysis_set, .data$exposure_category),
    label = factor(.data$label, levels = rev(paste(
      rep(analysis_set_order[-1], each = 2),
      rep(c("<35", ">50"), times = length(analysis_set_order[-1]))
    ))),
    hr_label = if_else(
      is.na(.data$adjusted_hr),
      "",
      paste0(sprintf("%.2f", .data$adjusted_hr), " (", sprintf("%.2f", .data$ci_low), "-", sprintf("%.2f", .data$ci_high), ")")
    )
  )

figure_s5 <- ggplot(figure_s5_data, aes(x = .data$adjusted_hr, y = .data$label, color = .data$analysis_set)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_errorbar(aes(xmin = .data$ci_low, xmax = .data$ci_high), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 2.75, label = .data$hr_label), hjust = 0, size = 2.8, color = "#222222") +
  scale_color_manual(values = figure_colors[analysis_set_order[-1]]) +
  scale_x_log10(limits = c(0.50, 4.4), breaks = c(0.75, 1, 1.5, 2, 3)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Adjusted hazard ratio vs 35-50 mmHg",
    y = NULL
  ) +
  base_theme(10) +
  theme(
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

# Result text ----------------------------------------------------------------

get_landmark_est <- function(analysis_set, category) {
  landmark_raw %>%
    filter(.data$analysis_set == .env$analysis_set, .data$exposure_category == .env$category) %>%
    slice(1)
}

get_burden_est <- function(metric, category) {
  table3_raw %>%
    filter(.data$exposure_metric == .env$metric, .data$category_or_contrast == .env$category) %>%
    slice(1)
}

pooled_lt35 <- get_landmark_est("Pooled", "<35")
pooled_gt50 <- get_landmark_est("Pooled", ">50")
mimic_lt35 <- get_landmark_est("MIMIC-IV", "<35")
mimic_gt50 <- get_landmark_est("MIMIC-IV", ">50")
amsterdam_lt35 <- get_landmark_est("AmsterdamUMCdb", "<35")
amsterdam_gt50 <- get_landmark_est("AmsterdamUMCdb", ">50")
chinese_lt35 <- get_landmark_est("Chinese", "<35")
chinese_gt50 <- get_landmark_est("Chinese", ">50")

any_yes <- get_burden_est("Any high-risk exposure", "Yes")
days_1 <- get_burden_est("Number of high-risk days", "1")
days_23 <- get_burden_est("Number of high-risk days", "2-3")
days_4 <- get_burden_est("Number of high-risk days", ">=4")
dose_per10 <- get_burden_est("Cumulative deviation dose", "per 10 mmHg-days")

any_exposed_patients <- ever_high_risk_summary$any_exposed_patients[1]
any_total_patients <- ever_high_risk_summary$n_patients[1]
any_exposed_pct <- ever_high_risk_summary$any_exposed_pct[1]

table_s6 <- tibble(
  `Population` = "Pooled cohort",
  `Support item` = "At least one high-risk PaCO2 exposure day during available ICU Days 1-7",
  `Denominator definition` = "Patients contributing at least one landmark day window",
  `Patients in denominator` = any_total_patients,
  `Patients meeting support item` = any_exposed_patients,
  `Percentage` = sprintf("%.1f%%", any_exposed_pct),
  `Derivation rule` = "For each global_stay_id, use the last observed landmark row in Result_3_high_risk_burden_landmark_qc.csv and count cumulative_high_risk_days > 0.",
  `Use limitation` = "Supports the descriptive burden sentence only; model estimates remain in Table 3."
)

get_fine_est <- function(analysis_set, category) {
  fine_raw %>%
    filter(.data$analysis_set == .env$analysis_set, .data$exposure_category == .env$category) %>%
    slice(1)
}

pooled_fine_lt35 <- get_fine_est("Pooled", "<35")
pooled_fine_35_40 <- get_fine_est("Pooled", "35-40")
pooled_fine_45_50 <- get_fine_est("Pooled", "45-50")
pooled_fine_gt50 <- get_fine_est("Pooled", ">50")

category_note_path <- file.path(MANUSCRIPT_DIR, "results_result3_category_sensitivity_note.txt")
category_note <- paste0(
  "PaCO2 五分类边界敏感性分析以 40-45 mmHg 作为参照。",
  "合并队列 stacked landmark 模型中，35-40 mmHg 和 45-50 mmHg 相比 40-45 mmHg 的 adjusted HR 分别为 ",
  fmt_estimate(pooled_fine_35_40$adjusted_hr, pooled_fine_35_40$ci_low, pooled_fine_35_40$ci_high),
  "和 ",
  fmt_estimate(pooled_fine_45_50$adjusted_hr, pooled_fine_45_50$ci_low, pooled_fine_45_50$ci_high),
  "；<35 mmHg 和 >50 mmHg 相比 40-45 mmHg 的 adjusted HR 分别为 ",
  fmt_estimate(pooled_fine_lt35$adjusted_hr, pooled_fine_lt35$ci_low, pooled_fine_lt35$ci_high),
  "和 ",
  fmt_estimate(pooled_fine_gt50$adjusted_hr, pooled_fine_gt50$ci_low, pooled_fine_gt50$ci_high),
  "（Table S3；Table S4）。"
)
writeLines(category_note, con = category_note_path, useBytes = TRUE)

result3_text <- c(
  "## Result 3. PaCO2 分类验证与高风险暴露负担",
  "",
  paste0(
    "以 35-50 mmHg 为参照的 stacked landmark 模型显示，合并队列中 TWA-PaCO2 <35 mmHg 和 >50 mmHg 均与较高的 28 天全因死亡风险相关，",
    "adjusted HR 分别为 ",
    fmt_estimate(pooled_lt35$adjusted_hr, pooled_lt35$ci_low, pooled_lt35$ci_high),
    "和 ",
    fmt_estimate(pooled_gt50$adjusted_hr, pooled_gt50$ci_low, pooled_gt50$ci_high),
    "。MIMIC-IV 和 AmsterdamUMCdb 队列的估计方向与合并结果一致；Chinese 队列的两项估计置信区间均跨越 1（Figure 2；Table 2；Figure S5）。"
  ),
  "",
  category_note,
  "",
  paste0(
    "高风险 PaCO2 暴露负担定义为截至各 landmark 日累计出现 TWA-PaCO2 <35 或 >50 mmHg 的日窗数及累计偏离剂量（Table S5）。",
    "在合并队列至少形成 1 个 landmark 日窗的 ",
    fmt_int(any_total_patients),
    " 例患者中，",
    fmt_int(any_exposed_patients),
    " 例（",
    sprintf("%.1f", any_exposed_pct),
    "%）在可观察的 ICU 首周内至少出现 1 天高风险 PaCO2 暴露。",
    "在 stacked landmark 累计暴露负担分析中，与尚无累计高风险暴露相比，任意累计高风险暴露的 adjusted HR 为 ",
    fmt_estimate(any_yes$adjusted_hr, any_yes$ci_low, any_yes$ci_high),
    "。按累计高风险暴露天数分层时，以 0 天为参照，1 天、2-3 天和 ≥4 天的 adjusted HR 分别为 ",
    fmt_estimate(days_1$adjusted_hr, days_1$ci_low, days_1$ci_high),
    "、",
    fmt_estimate(days_23$adjusted_hr, days_23$ci_low, days_23$ci_high),
    "和 ",
    fmt_estimate(days_4$adjusted_hr, days_4$ci_low, days_4$ci_high),
    "，趋势检验 ",
    format_p(days_trend_p),
    "。累计高风险偏离剂量每增加 10 mmHg-days 的 adjusted HR 为 ",
    fmt_estimate(dose_per10$adjusted_hr, dose_per10$ci_low, dose_per10$ci_high),
    "（Table 3）。"
  )
)

# Export ---------------------------------------------------------------------

write_csv(table2, file.path(TABLE_MAIN_DIR, "Table_2_primary_landmark_categories.csv"), na = "")
write.xlsx(table2, file.path(TABLE_MAIN_DIR, "Table_2_primary_landmark_categories.xlsx"), overwrite = TRUE)
write_csv(landmark_detail_raw, file.path(TABLE_SUPP_DIR, "Result_3_primary_landmark_categories_raw.csv"), na = "")
write.xlsx(landmark_detail_raw, file.path(TABLE_SUPP_DIR, "Result_3_primary_landmark_categories_raw.xlsx"), overwrite = TRUE)
write_csv(table_s3, file.path(TABLE_SUPP_DIR, "Table_S3_fine_paco2_category_sensitivity_landmark.csv"), na = "")
write.xlsx(table_s3, file.path(TABLE_SUPP_DIR, "Table_S3_fine_paco2_category_sensitivity_landmark.xlsx"), overwrite = TRUE)
write_csv(fine_detail_raw, file.path(TABLE_SUPP_DIR, "Result_3_fine_paco2_category_sensitivity_raw.csv"), na = "")
write.xlsx(fine_detail_raw, file.path(TABLE_SUPP_DIR, "Result_3_fine_paco2_category_sensitivity_raw.xlsx"), overwrite = TRUE)
write_csv(table_s4, file.path(TABLE_SUPP_DIR, "Table_S4_paco2_category_distribution_event_rates.csv"), na = "")
write.xlsx(table_s4, file.path(TABLE_SUPP_DIR, "Table_S4_paco2_category_distribution_event_rates.xlsx"), overwrite = TRUE)
write_csv(
  figure2_data,
  file.path(TABLE_SUPP_DIR, "Figure_2_day_specific_landmark_data.csv"),
  na = ""
)
write.xlsx(
  figure2_data,
  file.path(TABLE_SUPP_DIR, "Figure_2_day_specific_landmark_data.xlsx"),
  overwrite = TRUE
)
write_csv(
  figure_s5_data,
  file.path(TABLE_SUPP_DIR, "Figure_S5_cohort_day_specific_summary_data.csv"),
  na = ""
)
write.xlsx(
  figure_s5_data,
  file.path(TABLE_SUPP_DIR, "Figure_S5_cohort_day_specific_summary_data.xlsx"),
  overwrite = TRUE
)
write_csv(
  stacked_landmark_hierarchy_raw,
  file.path(TABLE_SUPP_DIR, "Result_3_stacked_pooled_landmark_adjustment_hierarchy_qc.csv"),
  na = ""
)
write_csv(
  stacked_cohort_hierarchy_raw,
  file.path(TABLE_SUPP_DIR, "Result_3_stacked_cohort_landmark_adjustment_hierarchy_qc.csv"),
  na = ""
)

write_csv(table3, file.path(TABLE_MAIN_DIR, "Table_3_high_risk_paco2_burden.csv"), na = "")
write.xlsx(table3, file.path(TABLE_MAIN_DIR, "Table_3_high_risk_paco2_burden.xlsx"), overwrite = TRUE)
write_csv(table3_raw, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_raw.csv"), na = "")
write.xlsx(table3_raw, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_raw.xlsx"), overwrite = TRUE)
write_csv(burden_day_raw, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_day_specific_raw.csv"), na = "")
write.xlsx(burden_day_raw, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_day_specific_raw.xlsx"), overwrite = TRUE)
write_csv(burden_counts_by_day, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_day_specific_counts.csv"), na = "")
write.xlsx(burden_counts_by_day, file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_day_specific_counts.xlsx"), overwrite = TRUE)
write_csv(
  burden_pooled %>%
    select(
      "global_stay_id",
      "analysis_cohort",
      "icu_day",
      "twa_paco2",
      "high_risk_day",
      "cumulative_high_risk_days",
      "cumulative_high_risk_prop",
      "cumulative_deviation_dose",
      "landmark_event",
      "landmark_time"
    ),
  file.path(TABLE_SUPP_DIR, "Result_3_high_risk_burden_landmark_qc.csv"),
  na = ""
)

write_csv(table_s5, file.path(TABLE_SUPP_DIR, "Table_S5_secondary_exposure_definitions.csv"), na = "")
write.xlsx(table_s5, file.path(TABLE_SUPP_DIR, "Table_S5_secondary_exposure_definitions.xlsx"), overwrite = TRUE)
write_csv(table_s6, file.path(TABLE_SUPP_DIR, "Table_S6_patient_level_high_risk_burden_support.csv"), na = "")
write.xlsx(table_s6, file.path(TABLE_SUPP_DIR, "Table_S6_patient_level_high_risk_burden_support.xlsx"), overwrite = TRUE)

ggsave(
  file.path(FIG_MAIN_DIR, "Figure_2_landmark_categorical_validation.png"),
  figure2,
  width = 7.8,
  height = 4.0,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_MAIN_DIR, "Figure_2_landmark_categorical_validation.pdf"),
  figure2,
  width = 7.8,
  height = 4.0,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S5_cohort_day_specific_landmark_summary.png"),
  figure_s5,
  width = 7.8,
  height = 4.8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S5_cohort_day_specific_landmark_summary.pdf"),
  figure_s5,
  width = 7.8,
  height = 4.8,
  bg = "white"
)

writeLines(
  result3_text,
  con = file.path(MANUSCRIPT_DIR, "simulated_results_result3_zh.md"),
  useBytes = TRUE
)

message("Done Result 3.")
message("Figure 2: ", file.path(FIG_MAIN_DIR, "Figure_2_landmark_categorical_validation.png"))
message("Figure S5: ", file.path(FIG_SUPP_DIR, "Figure_S5_cohort_day_specific_landmark_summary.png"))
message("Table 2: ", file.path(TABLE_MAIN_DIR, "Table_2_primary_landmark_categories.xlsx"))
message("Table S3: ", file.path(TABLE_SUPP_DIR, "Table_S3_fine_paco2_category_sensitivity_landmark.xlsx"))
message("Table S4: ", file.path(TABLE_SUPP_DIR, "Table_S4_paco2_category_distribution_event_rates.xlsx"))
message("Table 3: ", file.path(TABLE_MAIN_DIR, "Table_3_high_risk_paco2_burden.xlsx"))
message("Table S5: ", file.path(TABLE_SUPP_DIR, "Table_S5_secondary_exposure_definitions.xlsx"))
message("Table S6: ", file.path(TABLE_SUPP_DIR, "Table_S6_patient_level_high_risk_burden_support.xlsx"))
message("Result text: ", file.path(MANUSCRIPT_DIR, "simulated_results_result3_zh.md"))
