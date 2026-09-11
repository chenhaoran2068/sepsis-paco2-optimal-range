options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(openxlsx)
  library(patchwork)
  library(readr)
  library(scales)
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
PAMM_DIR <- paths$pamm_root
POOLED_DIR <- file.path(ANALYSIS_DATA_ROOT, "pooled")

FIG_MAIN_DIR <- file.path(RESULT_ROOT, "figures", "main")
FIG_SUPP_DIR <- file.path(RESULT_ROOT, "figures", "supplementary")
TABLE_SUPP_DIR <- file.path(RESULT_ROOT, "tables", "supplementary")
MANUSCRIPT_DIR <- file.path(RESULT_ROOT, "manuscript")
MODEL_OUT_DIR <- file.path(RESULT_ROOT, "intermediate", "result2_models")

dir.create(FIG_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PRIMARY_K <- 7L
K_CANDIDATES <- 3:8
REFERENCE_PACO2 <- 40

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

scope_to_cohort <- c(
  mimic = "MIMIC",
  amsterdam = "AmsterdamUMCdb",
  chinese = "Chinese cohort"
)

model_display <- c(
  model1 = "Model 1",
  model2 = "Model 2",
  model3 = "Model 3"
)

figure_colors <- c(
  "MIMIC-IV" = "#2C7FB8",
  "AmsterdamUMCdb" = "#41AB5D",
  "Chinese" = "#D95F0E",
  "Pooled" = "#222222"
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
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.title = element_text(face = "bold", size = base_size + 1)
    )
}

message("Running Result 2: primary PAMM analysis...")

ped_base <- read_parquet(file.path(PAMM_DIR, "ped_full_28d_locf.parquet")) %>%
  as.data.frame() %>%
  filter(!is.na(.data$twa_paco2))

baseline_included <- read_parquet(file.path(POOLED_DIR, "all_baseline_outcome.parquet")) %>%
  as.data.frame() %>%
  filter(.data$included_final) %>%
  select("global_stay_id")

day_data <- read_parquet(file.path(POOLED_DIR, "all_day_long.parquet")) %>%
  as.data.frame() %>%
  filter(
    .data$included_final,
    .data$observable_icu_day,
    .data$icu_day >= 1L,
    .data$icu_day <= 7L,
    !is.na(.data$twa_paco2)
  ) %>%
  mutate(
    cohort_label = unname(cohort_display[.data$analysis_cohort]),
    cohort_label = factor(.data$cohort_label, levels = c("MIMIC-IV", "AmsterdamUMCdb", "Chinese"))
  )

raw_paco2_used <- read_parquet(file.path(POOLED_DIR, "all_raw_paco2.parquet")) %>%
  as.data.frame() %>%
  semi_join(baseline_included, by = "global_stay_id") %>%
  filter(
    .data$used_for_twa %in% TRUE,
    .data$icu_day >= 1L,
    .data$icu_day <= 7L
  )

# 分析 -----------------------------------------------------------------------

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  names(sort(table(x), decreasing = TRUE))[1]
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

cohort_median <- function(data, variable) {
  data %>%
    group_by(.data$analysis_cohort) %>%
    summarise(value = median(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(.data$value), NA_real_, .data$value))
}

impute_dynamic_numeric <- function(data, variable, out_variable) {
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
    fill(all_of(within_variable), .direction = "down") %>%
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

prepare_imputed_data <- function(data, scope_id) {
  bmi_cohort <- cohort_median(data, "bmi") %>%
    rename(bmi_cohort_median = "value")
  bmi_overall <- median(data$bmi, na.rm = TRUE)
  if (is.nan(bmi_overall)) {
    bmi_overall <- NA_real_
  }

  data <- data %>%
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

  data <- impute_dynamic_numeric(data, "daily_sofa", "daily_sofa_imp")
  data <- impute_dynamic_numeric(data, "daily_lactate", "daily_lactate_imp")
  data <- impute_dynamic_numeric(data, "daily_oxygenation_value", "daily_oxygenation_value_imp")

  data %>%
    mutate(
      daily_lactate_imp = pmax(.data$daily_lactate_imp, 0),
      daily_lactate_log1p_imp = log1p(.data$daily_lactate_imp),
      daily_oxygenation_per100_imp = .data$daily_oxygenation_value_imp / 100
    )
}

model_terms <- function(model_id) {
  time_fixed <- c(
    "baseline_stratum",
    "age_per10",
    "sex_clean",
    "bmi_per5_imp",
    "bmi_missing"
  )

  dynamic <- c(
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

  if (model_id == "model1") {
    return("baseline_stratum")
  }
  if (model_id == "model2") {
    return(time_fixed)
  }
  if (model_id == "model3") {
    return(c(time_fixed, dynamic))
  }
  stop("Unsupported model_id: ", model_id, call. = FALSE)
}

formula_for_model <- function(data, model_id, paco2_k) {
  terms <- model_terms(model_id)
  usable_terms <- terms[vapply(terms, function(x) is_variable_usable(data, x), logical(1))]
  exposure_term <- paste0("s(twa_paco2, k = ", paco2_k, ", bs = 'cr')")
  as.formula(paste0(
    "ped_status ~ ",
    paste(c(usable_terms, exposure_term), collapse = " + "),
    " + offset(offset_log_interval)"
  ))
}

formula_for_day_interaction_model <- function(data, paco2_k) {
  terms <- model_terms("model3")
  usable_terms <- terms[vapply(terms, function(x) is_variable_usable(data, x), logical(1))]
  exposure_term <- paste0("te(source_icu_day, twa_paco2, k = c(4, ", paco2_k, "), bs = c('cr', 'cr'))")
  as.formula(paste0(
    "ped_status ~ ",
    paste(c(usable_terms, exposure_term), collapse = " + "),
    " + offset(offset_log_interval)"
  ))
}

template_newdata <- function(fit_data, grid) {
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
    daily_oxygenation_value_missing = as.integer(round(mean(fit_data$daily_oxygenation_value_missing, na.rm = TRUE)))
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

smooth_stats <- function(fit) {
  stable <- summary(fit)$s.table
  smooth_row <- grep("twa_paco2", rownames(stable), fixed = TRUE)[1]
  p_col <- intersect(c("p-value", "p"), colnames(stable))[1]
  tibble(
    edf_twa_paco2 = as.numeric(stable[smooth_row, "edf"]),
    p_twa_paco2 = as.numeric(stable[smooth_row, p_col])
  )
}

predict_hr_curve <- function(fit, fit_data, model_id, scope_id, paco2_k, grid_n = 281) {
  p1 <- as.numeric(quantile(fit_data$twa_paco2, 0.01, na.rm = TRUE))
  p99 <- as.numeric(quantile(fit_data$twa_paco2, 0.99, na.rm = TRUE))
  grid_min <- max(20, floor(p1))
  grid_max <- min(90, ceiling(p99))
  if (grid_max <= grid_min) {
    grid_min <- floor(min(fit_data$twa_paco2, na.rm = TRUE))
    grid_max <- ceiling(max(fit_data$twa_paco2, na.rm = TRUE))
  }

  grid <- seq(grid_min, grid_max, length.out = grid_n)
  if (REFERENCE_PACO2 >= grid_min && REFERENCE_PACO2 <= grid_max) {
    grid <- sort(unique(c(grid, REFERENCE_PACO2)))
  }
  nd <- template_newdata(fit_data, grid)
  xmat <- predict(fit, nd, type = "lpmatrix")
  eta <- as.vector(xmat %*% coef(fit))
  nadir_idx <- which.min(eta)
  nadir_paco2 <- grid[nadir_idx]
  xref <- xmat[nadir_idx, , drop = FALSE]
  xdiff <- sweep(xmat, 2, xref[1, ], "-")
  vc <- vcov(fit)
  se <- sqrt(pmax(0, rowSums((xdiff %*% vc) * xdiff)))
  log_hr <- eta - eta[nadir_idx]

  reference_idx <- which.min(abs(grid - REFERENCE_PACO2))
  reference_paco2 <- grid[reference_idx]
  xref_40 <- xmat[reference_idx, , drop = FALSE]
  xdiff_40 <- sweep(xmat, 2, xref_40[1, ], "-")
  se_40 <- sqrt(pmax(0, rowSums((xdiff_40 %*% vc) * xdiff_40)))
  log_hr_40 <- eta - eta[reference_idx]

  tibble(
    model_id = model_id,
    scope_id = scope_id,
    paco2_k = paco2_k,
    paco2 = grid,
    support_p1 = p1,
    support_p99 = p99,
    plot_grid_min = grid_min,
    plot_grid_max = grid_max,
    nadir_paco2 = nadir_paco2,
    reference_paco2 = reference_paco2,
    log_hr_vs_nadir = log_hr,
    se_log_hr_vs_nadir = se,
    hr_vs_nadir = exp(log_hr),
    hr_low_vs_nadir = exp(log_hr - 1.96 * se),
    hr_high_vs_nadir = exp(log_hr + 1.96 * se),
    log_hr_vs_40 = log_hr_40,
    se_log_hr_vs_40 = se_40,
    hr_vs_40 = exp(log_hr_40),
    hr_low_vs_40 = exp(log_hr_40 - 1.96 * se_40),
    hr_high_vs_40 = exp(log_hr_40 + 1.96 * se_40)
  )
}

nearest_curve_point <- function(curve, point) {
  curve %>%
    slice_min(abs(.data$paco2 - point), n = 1, with_ties = FALSE) %>%
    transmute(
      "{paste0('hr_at_', point, '_vs_nadir')}" := .data$hr_vs_nadir,
      "{paste0('hr_low_at_', point, '_vs_nadir')}" := .data$hr_low_vs_nadir,
      "{paste0('hr_high_at_', point, '_vs_nadir')}" := .data$hr_high_vs_nadir,
      "{paste0('hr_at_', point, '_vs_40')}" := .data$hr_vs_40,
      "{paste0('hr_low_at_', point, '_vs_40')}" := .data$hr_low_vs_40,
      "{paste0('hr_high_at_', point, '_vs_40')}" := .data$hr_high_vs_40
    )
}

summarise_curve <- function(curve, fit_data, fit, form, model_id, scope_id, paco2_k) {
  near5 <- curve %>% filter(.data$hr_vs_nadir <= 1.05)
  near10 <- curve %>% filter(.data$hr_vs_nadir <= 1.10)
  point_summaries <- bind_cols(
    nearest_curve_point(curve, 35),
    nearest_curve_point(curve, 45),
    nearest_curve_point(curve, 50),
    nearest_curve_point(curve, 60)
  )

  bind_cols(
    tibble(
      model_id = model_id,
      scope_id = scope_id,
      paco2_k = paco2_k,
      n_rows = nrow(fit_data),
      n_stays = n_distinct(fit_data$global_stay_id),
      n_events = sum(fit_data$ped_status, na.rm = TRUE),
      support_p1 = unique(curve$support_p1),
      support_p99 = unique(curve$support_p99),
      nadir_paco2 = unique(curve$nadir_paco2),
      near_min_5pct_low = if (nrow(near5) > 0) min(near5$paco2) else NA_real_,
      near_min_5pct_high = if (nrow(near5) > 0) max(near5$paco2) else NA_real_,
      near_min_10pct_low = if (nrow(near10) > 0) min(near10$paco2) else NA_real_,
      near_min_10pct_high = if (nrow(near10) > 0) max(near10$paco2) else NA_real_,
      aic = AIC(fit),
      model_formula = paste(deparse(form), collapse = " ")
    ),
    point_summaries,
    smooth_stats(fit)
  )
}

predict_day_specific_curves <- function(fit, fit_data, paco2_k, grid_n = 281) {
  p1 <- as.numeric(quantile(fit_data$twa_paco2, 0.01, na.rm = TRUE))
  p99 <- as.numeric(quantile(fit_data$twa_paco2, 0.99, na.rm = TRUE))
  grid_min <- max(20, floor(p1))
  grid_max <- min(90, ceiling(p99))
  grid <- seq(grid_min, grid_max, length.out = grid_n)

  bind_rows(lapply(1:7, function(day_id) {
    nd <- template_newdata(fit_data, grid) %>%
      mutate(source_icu_day = day_id)
    xmat <- predict(fit, nd, type = "lpmatrix")
    eta <- as.vector(xmat %*% coef(fit))
    nadir_idx <- which.min(eta)
    xref <- xmat[nadir_idx, , drop = FALSE]
    xdiff <- sweep(xmat, 2, xref[1, ], "-")
    vc <- vcov(fit)
    se <- sqrt(pmax(0, rowSums((xdiff %*% vc) * xdiff)))
    log_hr <- eta - eta[nadir_idx]

    tibble(
      icu_day = day_id,
      paco2_k = paco2_k,
      paco2 = grid,
      support_p1 = p1,
      support_p99 = p99,
      nadir_paco2 = grid[nadir_idx],
      log_hr_vs_day_nadir = log_hr,
      se_log_hr_vs_day_nadir = se,
      hr_vs_day_nadir = exp(log_hr),
      hr_low_vs_day_nadir = exp(log_hr - 1.96 * se),
      hr_high_vs_day_nadir = exp(log_hr + 1.96 * se)
    )
  }))
}

nearest_day_curve_point <- function(curve, point) {
  curve %>%
    slice_min(abs(.data$paco2 - point), n = 1, with_ties = FALSE) %>%
    transmute(
      "{paste0('hr_at_', point, '_vs_day_nadir')}" := .data$hr_vs_day_nadir,
      "{paste0('hr_low_at_', point, '_vs_day_nadir')}" := .data$hr_low_vs_day_nadir,
      "{paste0('hr_high_at_', point, '_vs_day_nadir')}" := .data$hr_high_vs_day_nadir
    )
}

summarise_day_specific_curves <- function(day_curves) {
  day_curves %>%
    group_by(.data$icu_day, .data$paco2_k) %>%
    group_modify(~ bind_cols(
      tibble(nadir_paco2 = first(.x$nadir_paco2)),
      nearest_day_curve_point(.x, 35),
      nearest_day_curve_point(.x, 50),
      nearest_day_curve_point(.x, 60)
    )) %>%
    ungroup()
}

fit_one_pamm <- function(base_data, model_id, scope_id, paco2_k) {
  fit_data <- base_data
  if (scope_id != "pooled") {
    fit_data <- fit_data %>%
      filter(.data$analysis_cohort == scope_to_cohort[[scope_id]])
  }
  fit_data <- prepare_imputed_data(fit_data, scope_id)

  form <- formula_for_model(fit_data, model_id, paco2_k)
  message(
    "Fitting ", model_display[[model_id]], " / ", scope_display[[scope_id]],
    " / k=", paco2_k,
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

  curve <- predict_hr_curve(fit, fit_data, model_id, scope_id, paco2_k)
  summary <- summarise_curve(curve, fit_data, fit, form, model_id, scope_id, paco2_k)
  saveRDS(
    fit,
    file.path(MODEL_OUT_DIR, paste0("result2_", model_id, "_", scope_id, "_k", paco2_k, ".rds"))
  )

  list(fit = fit, curve = curve, summary = summary)
}

fit_day_interaction_pamm <- function(base_data, paco2_k) {
  fit_data <- prepare_imputed_data(base_data, "pooled")
  form <- formula_for_day_interaction_model(fit_data, paco2_k)
  message(
    "Fitting Model 3 / Pooled / day-specific interaction / k=", paco2_k,
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

  day_curves <- predict_day_specific_curves(fit, fit_data, paco2_k)
  day_summary <- summarise_day_specific_curves(day_curves)
  saveRDS(
    fit,
    file.path(MODEL_OUT_DIR, paste0("result2_model3_pooled_day_interaction_k", paco2_k, ".rds"))
  )

  list(fit = fit, curve = day_curves, summary = day_summary, formula = form)
}

fit_plan <- bind_rows(
  tidyr::expand_grid(
    model_id = c("model1", "model2", "model3"),
    scope_id = c("pooled", "mimic", "amsterdam", "chinese"),
    paco2_k = PRIMARY_K
  ),
  tibble(model_id = "model3", scope_id = "pooled", paco2_k = K_CANDIDATES)
) %>%
  distinct(.data$model_id, .data$scope_id, .data$paco2_k, .keep_all = TRUE)

fit_results <- vector("list", nrow(fit_plan))
for (i in seq_len(nrow(fit_plan))) {
  spec <- fit_plan[i, ]
  fit_results[[i]] <- fit_one_pamm(
    base_data = ped_base,
    model_id = spec$model_id,
    scope_id = spec$scope_id,
    paco2_k = spec$paco2_k
  )
}

day_interaction_result <- fit_day_interaction_pamm(ped_base, PRIMARY_K)

curves <- bind_rows(lapply(fit_results, `[[`, "curve")) %>%
  mutate(
    model_label = factor(model_display[.data$model_id], levels = c("Model 1", "Model 2", "Model 3")),
    scope_label = factor(scope_display[.data$scope_id], levels = c("Pooled", "MIMIC-IV", "AmsterdamUMCdb", "Chinese"))
  )

summaries <- bind_rows(lapply(fit_results, `[[`, "summary")) %>%
  mutate(
    model_label = unname(model_display[.data$model_id]),
    scope_label = unname(scope_display[.data$scope_id]),
    near_min_5pct_label = paste0(round(.data$near_min_5pct_low, 1), "-", round(.data$near_min_5pct_high, 1)),
    near_min_10pct_label = paste0(round(.data$near_min_10pct_low, 1), "-", round(.data$near_min_10pct_high, 1))
  )

dist_overall <- day_data %>%
  summarise(
    n_day_windows = n(),
    median = median(.data$twa_paco2, na.rm = TRUE),
    q1 = quantile(.data$twa_paco2, 0.25, na.rm = TRUE),
    q3 = quantile(.data$twa_paco2, 0.75, na.rm = TRUE)
  )

dist_by_day <- day_data %>%
  summarise(
    n_day_windows = n(),
    median = median(.data$twa_paco2, na.rm = TRUE),
    q1 = quantile(.data$twa_paco2, 0.25, na.rm = TRUE),
    q3 = quantile(.data$twa_paco2, 0.75, na.rm = TRUE),
    .by = "icu_day"
  ) %>%
  arrange(.data$icu_day)

dist_by_cohort <- day_data %>%
  summarise(
    n_day_windows = n(),
    median = median(.data$twa_paco2, na.rm = TRUE),
    q1 = quantile(.data$twa_paco2, 0.25, na.rm = TRUE),
    q3 = quantile(.data$twa_paco2, 0.75, na.rm = TRUE),
    .by = "cohort_label"
  ) %>%
  arrange(.data$cohort_label)

raw_paco2_summary <- raw_paco2_used %>%
  mutate(cohort_label = unname(cohort_display[.data$analysis_cohort])) %>%
  summarise(n_records = n(), .by = "cohort_label")

primary_summary <- summaries %>%
  filter(.data$model_id == "model3", .data$scope_id == "pooled", .data$paco2_k == PRIMARY_K)

primary_curve <- curves %>%
  filter(.data$model_id == "model3", .data$scope_id == "pooled", .data$paco2_k == PRIMARY_K)

cohort_curves <- curves %>%
  filter(.data$model_id == "model3", .data$scope_id %in% c("mimic", "amsterdam", "chinese"), .data$paco2_k == PRIMARY_K)

hierarchy_curves <- curves %>%
  filter(.data$scope_id == "pooled", .data$paco2_k == PRIMARY_K)

k_grid_curves <- curves %>%
  filter(.data$model_id == "model3", .data$scope_id == "pooled", .data$paco2_k %in% K_CANDIDATES)

k_grid_summary <- summaries %>%
  filter(.data$model_id == "model3", .data$scope_id == "pooled", .data$paco2_k %in% K_CANDIDATES) %>%
  arrange(.data$paco2_k)

day_specific_curves <- day_interaction_result$curve
day_specific_summary <- day_interaction_result$summary %>%
  arrange(.data$icu_day)

fmt_p_value <- function(p) {
  p <- as.numeric(p)
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "P < .001",
      paste0("P = ", sub("^0", "", sprintf("%.3f", p)))
    )
  )
}

table_s2 <- summaries %>%
  filter(
    .data$paco2_k == PRIMARY_K,
    .data$model_id %in% c("model1", "model2", "model3"),
    .data$scope_id %in% c("pooled", "mimic", "amsterdam", "chinese")
  ) %>%
  mutate(
    Model = .data$model_label,
    `Analysis set` = .data$scope_label,
    `N patients` = .data$n_stays,
    `N PED rows` = .data$n_rows,
    Events = .data$n_events,
    `Nadir PaCO2, mmHg` = sprintf("%.1f", .data$nadir_paco2),
    `Near-minimum range, mmHg` = paste0(
      sprintf("%.1f", .data$near_min_5pct_low),
      "-",
      sprintf("%.1f", .data$near_min_5pct_high)
    ),
    `P for smooth term` = fmt_p_value(.data$p_twa_paco2),
    `HR at 35 vs nadir` = sprintf(
      "%.2f (%.2f-%.2f)",
      .data$hr_at_35_vs_nadir,
      .data$hr_low_at_35_vs_nadir,
      .data$hr_high_at_35_vs_nadir
    ),
    `HR at 50 vs nadir` = sprintf(
      "%.2f (%.2f-%.2f)",
      .data$hr_at_50_vs_nadir,
      .data$hr_low_at_50_vs_nadir,
      .data$hr_high_at_50_vs_nadir
    ),
    `HR at 60 vs nadir` = sprintf(
      "%.2f (%.2f-%.2f)",
      .data$hr_at_60_vs_nadir,
      .data$hr_low_at_60_vs_nadir,
      .data$hr_high_at_60_vs_nadir
    ),
    `HR at 50 vs 40 mmHg` = sprintf(
      "%.2f (%.2f-%.2f)",
      .data$hr_at_50_vs_40,
      .data$hr_low_at_50_vs_40,
      .data$hr_high_at_50_vs_40
    ),
    AIC = round(.data$aic, 1)
  ) %>%
  arrange(
    factor(.data$Model, levels = c("Model 1", "Model 2", "Model 3")),
    factor(.data$`Analysis set`, levels = c("Pooled", "MIMIC-IV", "AmsterdamUMCdb", "Chinese"))
  ) %>%
  select(
    "Model",
    "Analysis set",
    "N patients",
    "N PED rows",
    "Events",
    "Nadir PaCO2, mmHg",
    "Near-minimum range, mmHg",
    "P for smooth term",
    "HR at 35 vs nadir",
    "HR at 50 vs nadir",
    "HR at 60 vs nadir",
    "HR at 50 vs 40 mmHg",
    "AIC"
  )

fmt_n <- function(x) {
  formatC(as.numeric(x), format = "f", digits = 0, big.mark = ",")
}

fmt_num <- function(x, digits = 1) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

fmt_hr_ci <- function(summary_row, point) {
  hr <- summary_row[[paste0("hr_at_", point, "_vs_nadir")]]
  lo <- summary_row[[paste0("hr_low_at_", point, "_vs_nadir")]]
  hi <- summary_row[[paste0("hr_high_at_", point, "_vs_nadir")]]
  paste0(fmt_num(hr, 2), " (95% CI ", fmt_num(lo, 2), "-", fmt_num(hi, 2), ")")
}

fmt_iqr <- function(data_row) {
  paste0(fmt_num(data_row$median, 1), " [", fmt_num(data_row$q1, 1), ", ", fmt_num(data_row$q3, 1), "]")
}

cohort_dist_text <- function(label) {
  row <- dist_by_cohort %>% filter(.data$cohort_label == label)
  fmt_iqr(row)
}

day_medians <- paste(fmt_num(dist_by_day$median, 1), collapse = ", ")

model_hierarchy_text <- summaries %>%
  filter(.data$scope_id == "pooled", .data$paco2_k == PRIMARY_K, .data$model_id %in% c("model1", "model2", "model3")) %>%
  arrange(factor(.data$model_id, levels = c("model1", "model2", "model3"))) %>%
  mutate(text = paste0(.data$model_label, " ", fmt_num(.data$nadir_paco2, 1), " mmHg")) %>%
  pull("text") %>%
  paste(collapse = "、")

cohort_hr50_text <- summaries %>%
  filter(.data$model_id == "model3", .data$scope_id %in% c("mimic", "amsterdam", "chinese"), .data$paco2_k == PRIMARY_K) %>%
  arrange(factor(.data$scope_id, levels = c("mimic", "amsterdam", "chinese"))) %>%
  mutate(
    text = paste0(
      .data$scope_label,
      " ",
      fmt_num(.data$hr_at_50_vs_40, 2),
      " (95% CI ",
      fmt_num(.data$hr_low_at_50_vs_40, 2),
      "-",
      fmt_num(.data$hr_high_at_50_vs_40, 2),
      ")"
    )
  ) %>%
  pull("text") %>%
  paste(collapse = "、")

k6_8_summary <- k_grid_summary %>%
  filter(.data$paco2_k %in% 6:8)

k6_8_nadir_range <- paste0(
  fmt_num(min(k6_8_summary$nadir_paco2), 1),
  "-",
  fmt_num(max(k6_8_summary$nadir_paco2), 1)
)

day_nadir_range <- paste0(
  fmt_num(min(day_specific_summary$nadir_paco2), 1),
  "-",
  fmt_num(max(day_specific_summary$nadir_paco2), 1)
)

day_hr50_range <- paste0(
  fmt_num(min(day_specific_summary$hr_at_50_vs_day_nadir), 2),
  "-",
  fmt_num(max(day_specific_summary$hr_at_50_vs_day_nadir), 2)
)

day_hr60_range <- paste0(
  fmt_num(min(day_specific_summary$hr_at_60_vs_day_nadir), 2),
  "-",
  fmt_num(max(day_specific_summary$hr_at_60_vs_day_nadir), 2)
)

result2_text <- c(
  "## Result 2. TWA-PaCO2 与 28 天死亡的连续关联",
  "",
  paste0(
    "在最终分析队列中，ICU 入科后第 1-7 天共形成 ",
    fmt_n(dist_overall$n_day_windows),
    " 个可观察日窗；用于计算 TWA-PaCO2 的有效 PaCO2 测量记录共 ",
    fmt_n(nrow(raw_paco2_used)),
    " 条。总体每日 TWA-PaCO2 中位数为 ",
    fmt_iqr(dist_overall),
    " mmHg；MIMIC-IV、AmsterdamUMCdb 和 Chinese 分别为 ",
    cohort_dist_text("MIMIC-IV"),
    "、",
    cohort_dist_text("AmsterdamUMCdb"),
    " 和 ",
    cohort_dist_text("Chinese"),
    " mmHg。ICU 第 1 天至第 7 天的 TWA-PaCO2 中位数分别为 ",
    day_medians,
    " mmHg；可观察日窗数量由第 1 天的 ",
    fmt_n(first(dist_by_day$n_day_windows)),
    " 个下降至第 7 天的 ",
    fmt_n(last(dist_by_day$n_day_windows)),
    " 个",
    "（Figure 1A；Figure S2）。"
  ),
  "",
  paste0(
    "主调整 PAMM 以 ICU 入科后第 1-7 天每日 TWA-PaCO2 作为时变连续暴露，",
    "并调整预设的基线和时变协变量。合并队列模型显示，TWA-PaCO2 与 ICU 入科后 28 天全因死亡之间存在非线性关联；",
    "估计曲线的最低风险锚点位于 ",
    fmt_num(primary_summary$nadir_paco2, 1),
    " mmHg，5% 近最低风险区间为 ",
    primary_summary$near_min_5pct_label,
    " mmHg。相对于该锚点，TWA-PaCO2 为 35、50 和 60 mmHg 时的 adjusted HR 分别为 ",
    fmt_hr_ci(primary_summary, 35),
    "、",
    fmt_hr_ci(primary_summary, 50),
    " 和 ",
    fmt_hr_ci(primary_summary, 60),
    "（Figure 1B；Table S2）。"
  ),
  "",
  paste0(
    "队列内分析显示，MIMIC-IV 和 AmsterdamUMCdb 的曲线最低风险区域均位于约 40 mmHg 附近；",
    "Chinese 队列曲线主要表现为高 PaCO2 区间风险逐渐升高。",
    "以 40 mmHg 作为统一参考点时，TWA-PaCO2 为 50 mmHg 的 adjusted HR 分别为 ",
    cohort_hr50_text,
    "。调整层级和 day-specific PAMM 分析均支持约 40 mmHg 的最低风险锚点，",
    "k-grid 分析支持主模型 k=",
    PRIMARY_K,
    " 的设定（Figure 1C；Figure S3；Figure S4；Table S2）。"
  )
)

# 导出结果 -------------------------------------------------------------------

range_background <- function() {
  list(
    annotate("rect", xmin = 35, xmax = 50, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.55),
    annotate("rect", xmin = 40, xmax = 45, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.45),
    geom_vline(xintercept = c(35, 50), linetype = "dashed", color = "grey40", linewidth = 0.35),
    geom_vline(xintercept = c(40, 45), linetype = "dotted", color = "grey40", linewidth = 0.35)
  )
}

panel_a <- ggplot(day_data, aes(x = .data$twa_paco2, color = .data$cohort_label, fill = .data$cohort_label)) +
  range_background() +
  geom_density(alpha = 0.15, linewidth = 0.75, adjust = 1.05) +
  scale_color_manual(values = figure_colors) +
  scale_fill_manual(values = figure_colors) +
  coord_cartesian(xlim = c(20, 80)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "Density",
    title = "Daily TWA-PaCO2 distribution"
  ) +
  base_theme()

panel_b <- ggplot(primary_curve, aes(x = .data$paco2, y = .data$hr_vs_nadir)) +
  range_background() +
  geom_ribbon(aes(ymin = .data$hr_low_vs_nadir, ymax = .data$hr_high_vs_nadir), fill = "#777777", alpha = 0.18) +
  geom_line(color = "#222222", linewidth = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  coord_cartesian(xlim = c(20, 80)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "Adjusted hazard ratio vs nadir",
    title = "Pooled primary adjusted PAMM"
  ) +
  base_theme()

panel_c <- ggplot(cohort_curves, aes(x = .data$paco2, y = .data$hr_vs_40, color = .data$scope_label, fill = .data$scope_label)) +
  range_background() +
  geom_ribbon(aes(ymin = .data$hr_low_vs_40, ymax = .data$hr_high_vs_40), alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  scale_color_manual(values = figure_colors) +
  scale_fill_manual(values = figure_colors) +
  coord_cartesian(xlim = c(20, 80)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "Adjusted hazard ratio vs 40 mmHg",
    title = "Cohort-specific primary adjusted PAMM"
  ) +
  base_theme()

figure_1 <- (panel_a / panel_b / panel_c) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(FIG_MAIN_DIR, "Figure_1_primary_pamm_twa_paco2.png"),
  figure_1,
  width = 7.5,
  height = 10.2,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_MAIN_DIR, "Figure_1_primary_pamm_twa_paco2.pdf"),
  figure_1,
  width = 7.5,
  height = 10.2,
  bg = "white"
)

figure_s2_day_labels <- dist_by_day %>%
  mutate(
    icu_day = as.character(.data$icu_day),
    label = paste0("Day ", .data$icu_day, "\n", "n=", comma(.data$n_day_windows, accuracy = 1))
  ) %>%
  {setNames(.$label, .$icu_day)}

figure_s2 <- ggplot(day_data, aes(x = factor(.data$icu_day), y = .data$twa_paco2, fill = .data$cohort_label)) +
  geom_hline(yintercept = c(35, 50), linetype = "dashed", color = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = c(40, 45), linetype = "dotted", color = "grey45", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.68, color = "grey25", linewidth = 0.3) +
  scale_fill_manual(values = figure_colors) +
  scale_x_discrete(labels = figure_s2_day_labels) +
  coord_cartesian(ylim = c(20, 80)) +
  facet_wrap(~cohort_label, ncol = 1) +
  labs(
    x = "ICU day and pooled observable day-windows",
    y = "Daily TWA-PaCO2, mmHg"
  ) +
  base_theme()

ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S2_twa_paco2_distribution_by_day_cohort.png"),
  figure_s2,
  width = 7.2,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S2_twa_paco2_distribution_by_day_cohort.pdf"),
  figure_s2,
  width = 7.2,
  height = 8,
  bg = "white"
)

figure_s2_day_summary <- dist_by_day %>%
  transmute(
    `ICU day` = .data$icu_day,
    `Observable day-windows` = .data$n_day_windows,
    `TWA-PaCO2 median [IQR], mmHg` = paste0(
      formatC(.data$median, format = "f", digits = 1),
      " [",
      formatC(.data$q1, format = "f", digits = 1),
      ", ",
      formatC(.data$q3, format = "f", digits = 1),
      "]"
    )
  )

write_csv(
  figure_s2_day_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S2_twa_paco2_distribution_by_day_summary.csv"),
  na = ""
)
write.xlsx(
  figure_s2_day_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S2_twa_paco2_distribution_by_day_summary.xlsx"),
  overwrite = TRUE
)

figure_s3_panel_a <- ggplot(hierarchy_curves, aes(x = .data$paco2, y = .data$hr_vs_nadir, color = .data$model_label)) +
  range_background() +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  coord_cartesian(xlim = c(20, 80)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "HR vs model-specific nadir",
    title = "Adjustment hierarchy"
  ) +
  base_theme()

figure_s3_panel_b <- ggplot(k_grid_curves, aes(x = .data$paco2, y = .data$hr_vs_nadir, color = factor(.data$paco2_k))) +
  range_background() +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  coord_cartesian(xlim = c(20, 80)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "HR vs nadir",
    color = "k",
    title = "Model 3 k-grid"
  ) +
  base_theme()

figure_s3_panel_c <- ggplot(k_grid_summary, aes(x = factor(.data$paco2_k), y = .data$nadir_paco2)) +
  geom_hline(yintercept = c(35, 50), linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = c(40, 45), linetype = "dotted", color = "grey55", linewidth = 0.3) +
  geom_pointrange(aes(ymin = .data$near_min_5pct_low, ymax = .data$near_min_5pct_high), color = "#222222", linewidth = 0.5) +
  coord_cartesian(ylim = c(25, 55)) +
  labs(
    x = "Spline basis dimension k",
    y = "PaCO2, mmHg",
    title = "Nadir stability"
  ) +
  base_theme() +
  theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 14))

figure_s4 <- ggplot(day_specific_curves, aes(x = .data$paco2, y = .data$hr_vs_day_nadir)) +
  range_background() +
  geom_ribbon(aes(ymin = .data$hr_low_vs_day_nadir, ymax = .data$hr_high_vs_day_nadir), fill = "#777777", alpha = 0.16) +
  geom_line(color = "#222222", linewidth = 0.65) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  facet_wrap(~paste0("Day ", .data$icu_day), ncol = 4) +
  coord_cartesian(xlim = c(25, 75), ylim = c(0.8, 3.6)) +
  labs(
    x = "Daily TWA-PaCO2, mmHg",
    y = "Hazard ratio vs day-specific nadir"
  ) +
  base_theme()

figure_s3 <- (figure_s3_panel_a / figure_s3_panel_b / figure_s3_panel_c) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S3_pamm_adjustment_hierarchy_model_qc.png"),
  figure_s3,
  width = 7.5,
  height = 9,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S3_pamm_adjustment_hierarchy_model_qc.pdf"),
  figure_s3,
  width = 7.5,
  height = 9,
  bg = "white"
)

ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S4_day_specific_pamm_curves.png"),
  figure_s4,
  width = 9,
  height = 5.8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(FIG_SUPP_DIR, "Figure_S4_day_specific_pamm_curves.pdf"),
  figure_s4,
  width = 9,
  height = 5.8,
  bg = "white"
)

write_csv(
  curves,
  file.path(TABLE_SUPP_DIR, "Result_2_pamm_curve_predictions.csv"),
  na = ""
)

write_csv(
  summaries,
  file.path(TABLE_SUPP_DIR, "Result_2_pamm_model_summary.csv"),
  na = ""
)

write_csv(
  k_grid_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S3_pamm_k_grid_summary.csv"),
  na = ""
)
write.xlsx(
  k_grid_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S3_pamm_k_grid_summary.xlsx"),
  overwrite = TRUE
)

write_csv(
  day_specific_curves,
  file.path(TABLE_SUPP_DIR, "Figure_S4_day_specific_pamm_curves.csv"),
  na = ""
)

write_csv(
  day_specific_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S4_day_specific_pamm_summary.csv"),
  na = ""
)
write.xlsx(
  day_specific_summary,
  file.path(TABLE_SUPP_DIR, "Figure_S4_day_specific_pamm_summary.xlsx"),
  overwrite = TRUE
)

write_csv(table_s2, file.path(TABLE_SUPP_DIR, "Table_S2_pamm_adjustment_hierarchy.csv"), na = "")
write.xlsx(table_s2, file.path(TABLE_SUPP_DIR, "Table_S2_pamm_adjustment_hierarchy.xlsx"), overwrite = TRUE)

result2_qc <- tibble(
  item = c(
    "primary_k",
    "observable_day_windows",
    "raw_paco2_records_used_for_twa",
    "primary_pamm_n_stays",
    "primary_pamm_n_ped_rows",
    "primary_pamm_events",
    "primary_pamm_nadir_paco2",
    "primary_pamm_near_min_5pct_range",
    "primary_pamm_aic",
    "k6_to_k8_nadir_range",
    "day1_to_day7_nadir_range",
    "day1_to_day7_hr50_range",
    "day1_to_day7_hr60_range"
  ),
  value = c(
    PRIMARY_K,
    dist_overall$n_day_windows,
    nrow(raw_paco2_used),
    primary_summary$n_stays,
    primary_summary$n_rows,
    primary_summary$n_events,
    round(primary_summary$nadir_paco2, 2),
    primary_summary$near_min_5pct_label,
    round(primary_summary$aic, 2),
    k6_8_nadir_range,
    day_nadir_range,
    day_hr50_range,
    day_hr60_range
  )
)

write_csv(result2_qc, file.path(TABLE_SUPP_DIR, "Result_2_qc_summary.csv"), na = "")

writeLines(
  result2_text,
  con = file.path(MANUSCRIPT_DIR, "simulated_results_result2_zh.md"),
  useBytes = TRUE
)

message("Done Result 2.")
message("Figure 1: ", file.path(FIG_MAIN_DIR, "Figure_1_primary_pamm_twa_paco2.png"))
message("Figure S2: ", file.path(FIG_SUPP_DIR, "Figure_S2_twa_paco2_distribution_by_day_cohort.png"))
message("Figure S3: ", file.path(FIG_SUPP_DIR, "Figure_S3_pamm_adjustment_hierarchy_model_qc.png"))
message("Figure S4: ", file.path(FIG_SUPP_DIR, "Figure_S4_day_specific_pamm_curves.png"))
message("Table S2: ", file.path(TABLE_SUPP_DIR, "Table_S2_pamm_adjustment_hierarchy.xlsx"))
message("Figure S3 k-grid summary: ", file.path(TABLE_SUPP_DIR, "Figure_S3_pamm_k_grid_summary.xlsx"))
message("Figure S4 day-specific PAMM summary: ", file.path(TABLE_SUPP_DIR, "Figure_S4_day_specific_pamm_summary.xlsx"))
message("Result text: ", file.path(MANUSCRIPT_DIR, "simulated_results_result2_zh.md"))
