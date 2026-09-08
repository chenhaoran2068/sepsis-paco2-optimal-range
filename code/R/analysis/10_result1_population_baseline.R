options(stringsAsFactors = FALSE)

# 配置环境（载入包、载入数据） -----------------------------------------------

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(openxlsx)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
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
DATA_DIR <- file.path(ANALYSIS_DATA_ROOT, "pooled")
TEMPLATE_DIR <- file.path(SCRIPT_ROOT, "assets", "templates")
TABLE_MAIN_DIR <- file.path(RESULT_ROOT, "tables", "main")
TABLE_SUPP_DIR <- file.path(RESULT_ROOT, "tables", "supplementary")
FIG_SUPP_DIR <- file.path(RESULT_ROOT, "figures", "supplementary")
MANUSCRIPT_DIR <- file.path(RESULT_ROOT, "manuscript")

dir.create(TABLE_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Running Result 1: population, baseline characteristics, and missingness...")

baseline <- read_parquet(file.path(DATA_DIR, "all_baseline_outcome.parquet")) %>%
  as.data.frame()

day_long <- read_parquet(file.path(DATA_DIR, "all_day_long.parquet")) %>%
  as.data.frame()

table1_template <- read_csv(
  file.path(TEMPLATE_DIR, "Table_1_baseline_characteristics_template.csv"),
  show_col_types = FALSE
)

table_s1_template <- read_csv(
  file.path(TEMPLATE_DIR, "Table_S1_missingness_template.csv"),
  show_col_types = FALSE
)

cohort_levels <- c("MIMIC", "AmsterdamUMCdb", "Chinese cohort")
cohort_display <- c(
  MIMIC = "MIMIC-IV",
  AmsterdamUMCdb = "AmsterdamUMCdb",
  `Chinese cohort` = "Chinese"
)

included_baseline <- baseline %>%
  filter(.data$included_final) %>%
  mutate(
    analysis_cohort = factor(.data$analysis_cohort, levels = cohort_levels),
    sex_clean = case_when(
      tolower(.data$sex) %in% c("male", "m") ~ "male",
      tolower(.data$sex) %in% c("female", "f") ~ "female",
      TRUE ~ NA_character_
    )
  )

included_day <- day_long %>%
  filter(.data$included_final, .data$observable_icu_day) %>%
  mutate(analysis_cohort = factor(.data$analysis_cohort, levels = cohort_levels))

included_day1 <- included_day %>%
  filter(.data$icu_day == 1L)

# 分析 -----------------------------------------------------------------------

fmt_n <- function(x) {
  comma(as.numeric(x), accuracy = 1)
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "", paste0(formatC(x, format = "f", digits = digits), "%"))
}

fmt_n_pct <- function(n, denom, digits = 1) {
  if (is.na(denom) || denom == 0) {
    return("")
  }
  paste0(fmt_n(n), " (", fmt_pct(n / denom * 100, digits), ")")
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return("")
  }
  paste0(
    formatC(mean(x), format = "f", digits = digits),
    " (",
    formatC(sd(x), format = "f", digits = digits),
    ")"
  )
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return("")
  }
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 2)
  paste0(
    formatC(q[2], format = "f", digits = digits),
    " [",
    formatC(q[1], format = "f", digits = digits),
    ", ",
    formatC(q[3], format = "f", digits = digits),
    "]"
  )
}

continuous_value <- function(data, variable, display, digits = 1) {
  if (display == "mean (SD)") {
    return(fmt_mean_sd(data[[variable]], digits))
  }
  if (display == "median [IQR]") {
    return(fmt_median_iqr(data[[variable]], digits))
  }
  stop("Unsupported continuous display: ", display, call. = FALSE)
}

categorical_value <- function(data, variable, positive_value = TRUE) {
  x <- data[[variable]]
  denom <- sum(!is.na(x))
  n <- sum(x %in% positive_value, na.rm = TRUE)
  fmt_n_pct(n, denom)
}

male_value <- function(data) {
  denom <- sum(!is.na(data$sex_clean))
  n <- sum(data$sex_clean == "male", na.rm = TRUE)
  fmt_n_pct(n, denom)
}

row_value <- function(data, characteristic, display) {
  switch(
    characteristic,
    "Study population" = fmt_n(nrow(data)),
    "28-day all-cause mortality" = categorical_value(data, "death_28d"),
    "Age" = continuous_value(data, "age", display, digits = 1),
    "Male sex" = male_value(data),
    "BMI" = continuous_value(data, "bmi", display, digits = 1),
    "COPD" = categorical_value(data, "copd"),
    "Baseline SOFA" = continuous_value(data, "baseline_sofa", display, digits = 1),
    "Baseline lactate" = continuous_value(data, "baseline_lactate", display, digits = 1),
    "Baseline pH" = continuous_value(data, "baseline_ph", display, digits = 2),
    "Baseline HCO3" = continuous_value(data, "baseline_hco3", display, digits = 1),
    "Baseline PaO2/FiO2" = continuous_value(data, "baseline_pf_ratio", display, digits = 1),
    "Mechanical ventilation within 24 h" = categorical_value(data, "baseline_mv"),
    "Vasopressor use within 24 h" = categorical_value(data, "baseline_vasopressor"),
    "CRRT/RRT within 24 h" = categorical_value(data, "baseline_crrt_rrt"),
    stop("Unsupported Table 1 characteristic: ", characteristic, call. = FALSE)
  )
}

day1_row_value <- function(data, characteristic, display) {
  switch(
    characteristic,
    "Day 1 TWA-PaCO2" = continuous_value(data, "twa_paco2", display, digits = 1),
    stop("Unsupported day 1 characteristic: ", characteristic, call. = FALSE)
  )
}

table1 <- table1_template
for (i in seq_len(nrow(table1))) {
  characteristic <- table1$Characteristic[i]
  display <- table1$Display[i]

  if (characteristic == "Day 1 TWA-PaCO2") {
    source_data <- included_day1
    get_value <- day1_row_value
  } else {
    source_data <- included_baseline
    get_value <- row_value
  }

  table1$Overall[i] <- get_value(source_data, characteristic, display)
  for (cohort in cohort_levels) {
    col_name <- cohort_display[[cohort]]
    table1[[col_name]][i] <- get_value(
      source_data %>% filter(.data$analysis_cohort == cohort),
      characteristic,
      display
    )
  }
}

table1_display_labels <- c(
  "Age" = "Age, years",
  "BMI" = "BMI, kg/m2",
  "Baseline lactate" = "Baseline lactate, mmol/L",
  "Day 1 TWA-PaCO2" = "Day 1 TWA-PaCO2, mmHg",
  "Baseline HCO3" = "Baseline HCO3, mmol/L",
  "Baseline PaO2/FiO2" = "Baseline PaO2/FiO2, mmHg"
)

table1 <- table1 %>%
  mutate(
    Characteristic = if_else(
      .data$Characteristic %in% names(table1_display_labels),
      unname(table1_display_labels[.data$Characteristic]),
      .data$Characteristic
    )
  ) %>%
  select(-any_of("Notes"))

missing_format <- function(data, variable) {
  denom <- nrow(data)
  n_missing <- sum(is.na(data[[variable]]))
  pct_missing <- n_missing / denom * 100
  if (n_missing > 0 && pct_missing < 0.1) {
    return(paste0(fmt_n(n_missing), " (<0.1%)"))
  }
  fmt_n_pct(n_missing, denom)
}

missing_by_cohort <- function(data, variable) {
  out <- tibble(
    Overall = missing_format(data, variable)
  )
  for (cohort in cohort_levels) {
    out[[cohort_display[[cohort]]]] <- missing_format(
      data %>% filter(.data$analysis_cohort == cohort),
      variable
    )
  }
  out
}

missing_variable_map <- list(
  "Age" = list(data = "baseline", variable = "age"),
  "Sex" = list(data = "baseline", variable = "sex"),
  "BMI" = list(data = "baseline", variable = "bmi"),
  "COPD" = list(data = "baseline", variable = "copd"),
  "Baseline SOFA" = list(data = "baseline", variable = "baseline_sofa"),
  "Daily SOFA days 1-7" = list(data = "day", variable = "daily_sofa"),
  "Baseline lactate" = list(data = "baseline", variable = "baseline_lactate"),
  "Daily lactate days 1-7" = list(data = "day", variable = "daily_lactate"),
  "Day 1 TWA-PaCO2" = list(data = "day1", variable = "twa_paco2"),
  "Daily TWA-PaCO2 days 1-7" = list(data = "day", variable = "twa_paco2"),
  "Baseline pH" = list(data = "baseline", variable = "baseline_ph"),
  "Daily pH days 1-7" = list(data = "day", variable = "daily_ph"),
  "Baseline HCO3" = list(data = "baseline", variable = "baseline_hco3"),
  "Baseline PaO2/FiO2" = list(data = "baseline", variable = "baseline_pf_ratio"),
  "Daily oxygenation measure days 1-7" = list(data = "day", variable = "daily_oxygenation_value"),
  "Mechanical ventilation within 24 h" = list(data = "baseline", variable = "baseline_mv"),
  "Daily mechanical ventilation days 1-7" = list(data = "day", variable = "daily_mv"),
  "Vasopressor use within 24 h" = list(data = "baseline", variable = "baseline_vasopressor"),
  "Daily vasopressor use days 1-7" = list(data = "day", variable = "daily_vasopressor"),
  "CRRT/RRT within 24 h" = list(data = "baseline", variable = "baseline_crrt_rrt"),
  "Daily CRRT/RRT days 1-7" = list(data = "day", variable = "daily_crrt_rrt"),
  "28-day all-cause mortality" = list(data = "baseline", variable = "death_28d")
)

get_missing_data <- function(data_name) {
  switch(
    data_name,
    baseline = included_baseline,
    day = included_day,
    day1 = included_day1,
    stop("Unsupported missingness data source: ", data_name, call. = FALSE)
  )
}

table_s1 <- table_s1_template
for (i in seq_len(nrow(table_s1))) {
  variable_label <- table_s1$Variable[i]
  spec <- missing_variable_map[[variable_label]]
  if (is.null(spec)) {
    stop("Missingness variable not mapped: ", variable_label, call. = FALSE)
  }
  values <- missing_by_cohort(get_missing_data(spec$data), spec$variable)
  table_s1$`Overall missing n (%)`[i] <- values$Overall
  table_s1$`MIMIC-IV missing n (%)`[i] <- values$`MIMIC-IV`
  table_s1$`AmsterdamUMCdb missing n (%)`[i] <- values$AmsterdamUMCdb
  table_s1$`Chinese missing n (%)`[i] <- values$Chinese
}

table_s1 <- table_s1 %>%
  select(
    "Variable",
    "Overall missing n (%)",
    "MIMIC-IV missing n (%)",
    "AmsterdamUMCdb missing n (%)",
    "Chinese missing n (%)"
  )

selection_steps <- tibble(
  step_order = 1:9,
  step_id = c(
    "source_icu_stays",
    "age_ge_18",
    "sepsis3_onset_window",
    "first_hospital_admission",
    "first_icu_admission",
    "icu_los_ge_24h",
    "outcome_28d_available",
    "paco2_complete_for_all_observable_icu_days",
    "final_analysis_cohort"
  ),
  label = c(
    "Raw ICU stays",
    "Age >=18 years",
    "Sepsis-3 onset\n-6 h to +24 h",
    "First hospitalization",
    "First ICU admission",
    "ICU LOS >=24 h",
    "28-day outcome\navailable",
    "All observable ICU days\nwith PaCO2",
    "Final cohort"
  )
)

cohort_flow_one <- function(data, cohort) {
  x <- data %>% filter(.data$analysis_cohort == cohort)
  remaining <- nrow(x)
  rows <- list(
    tibble(
      analysis_cohort = cohort,
      step_order = 1L,
      step_id = "source_icu_stays",
      n_before_step = NA_integer_,
      n_excluded_this_step = NA_integer_,
      n_remaining_after_step = remaining
    )
  )

  filters <- list(
    age_ge_18 = quote(.data$adult_flag),
    sepsis3_onset_window = quote(.data$sepsis3_window_flag),
    first_hospital_admission = quote(.data$first_hosp_flag),
    first_icu_admission = quote(.data$first_icu_flag),
    icu_los_ge_24h = quote(.data$icu_los_ge_24h),
    outcome_28d_available = quote(.data$outcome_28d_available),
    paco2_complete_for_all_observable_icu_days = quote(!.data$missing_observable_day_paco2_flag),
    final_analysis_cohort = quote(.data$included_final)
  )

  current <- x
  for (j in seq_along(filters)) {
    step_id <- names(filters)[j]
    before <- nrow(current)
    current <- current %>% filter(!!filters[[j]])
    after <- nrow(current)
    rows[[j + 1L]] <- tibble(
      analysis_cohort = cohort,
      step_order = j + 1L,
      step_id = step_id,
      n_before_step = before,
      n_excluded_this_step = before - after,
      n_remaining_after_step = after
    )
  }

  bind_rows(rows)
}

cohort_flow <- bind_rows(lapply(cohort_levels, function(x) cohort_flow_one(baseline, x))) %>%
  left_join(selection_steps, by = c("step_order", "step_id")) %>%
  mutate(
    cohort_label = recode(.data$analysis_cohort, !!!cohort_display),
    n_remaining_label = paste0("n = ", fmt_n(.data$n_remaining_after_step)),
    excluded_label = if_else(
      .data$step_order == 1L,
      "",
      paste0("Excluded: ", fmt_n(.data$n_excluded_this_step))
    )
  )

manuscript_stats <- function() {
  all_data <- included_baseline
  day1_all <- included_day1

  stat_cohort_n <- included_baseline %>%
    summarise(n = n(), deaths = sum(.data$death_28d %in% TRUE), .by = "analysis_cohort") %>%
    mutate(mortality = .data$deaths / .data$n * 100)

  get_n <- function(cohort) {
    stat_cohort_n %>% filter(.data$analysis_cohort == cohort) %>% pull(.data$n)
  }
  get_deaths <- function(cohort) {
    stat_cohort_n %>% filter(.data$analysis_cohort == cohort) %>% pull(.data$deaths)
  }
  get_mortality <- function(cohort) {
    stat_cohort_n %>% filter(.data$analysis_cohort == cohort) %>% pull(.data$mortality)
  }
  get_median <- function(data, var, digits = 1) {
    fmt_median_iqr(data[[var]], digits)
  }
  get_pct <- function(data, var) {
    categorical_value(data, var)
  }
  get_male_pct <- function(data) {
    male_value(data)
  }

  text <- c(
    "## Result 1. 研究人群与基线特征",
    "",
    paste0(
      "按纳入排除标准，本研究在 MIMIC-IV、AmsterdamUMCdb 和 Chinese 三个队列中识别 ICU 脓毒症患者，最终纳入 ",
      fmt_n(nrow(all_data)),
      " 例患者，其中 MIMIC-IV、AmsterdamUMCdb 和 Chinese 分别为 ",
      fmt_n(get_n("MIMIC")),
      "、",
      fmt_n(get_n("AmsterdamUMCdb")),
      " 和 ",
      fmt_n(get_n("Chinese cohort")),
      " 例（Figure S1）。"
    ),
    "",
    paste0(
      "总体 28 天全因死亡为 ",
      fmt_n(sum(all_data$death_28d %in% TRUE)),
      " 例（",
      fmt_pct(sum(all_data$death_28d %in% TRUE) / nrow(all_data) * 100),
      "）；年龄为 ",
      get_median(all_data, "age"),
      " 岁，男性为 ",
      get_male_pct(all_data),
      "，BMI 为 ",
      get_median(all_data, "bmi"),
      " kg/m2。MIMIC-IV、AmsterdamUMCdb 和 Chinese 的 28 天全因死亡分别为 ",
      fmt_n(get_deaths("MIMIC")),
      " 例（",
      fmt_pct(get_mortality("MIMIC")),
      "）、",
      fmt_n(get_deaths("AmsterdamUMCdb")),
      " 例（",
      fmt_pct(get_mortality("AmsterdamUMCdb")),
      "）和 ",
      fmt_n(get_deaths("Chinese cohort")),
      " 例（",
      fmt_pct(get_mortality("Chinese cohort")),
      "），年龄分别为 ",
      get_median(all_data %>% filter(.data$analysis_cohort == "MIMIC"), "age"),
      "、",
      get_median(all_data %>% filter(.data$analysis_cohort == "AmsterdamUMCdb"), "age"),
      " 和 ",
      get_median(all_data %>% filter(.data$analysis_cohort == "Chinese cohort"), "age"),
      " 岁，男性比例分别为 ",
      get_male_pct(all_data %>% filter(.data$analysis_cohort == "MIMIC")),
      "、",
      get_male_pct(all_data %>% filter(.data$analysis_cohort == "AmsterdamUMCdb")),
      " 和 ",
      get_male_pct(all_data %>% filter(.data$analysis_cohort == "Chinese cohort")),
      "。总体基线 SOFA 评分、乳酸、Day 1 TWA-PaCO2、pH、HCO3 和 PaO2/FiO2 分别为 ",
      get_median(all_data, "baseline_sofa"),
      "、",
      get_median(all_data, "baseline_lactate"),
      " mmol/L、",
      get_median(day1_all, "twa_paco2"),
      " mmHg、",
      get_median(all_data, "baseline_ph", digits = 2),
      "、",
      get_median(all_data, "baseline_hco3"),
      " mmol/L 和 ",
      get_median(all_data, "baseline_pf_ratio"),
      " mmHg。ICU 入科早期接受机械通气、血管活性药物和 CRRT/RRT 的比例分别为 ",
      get_pct(all_data, "baseline_mv"),
      "、",
      get_pct(all_data, "baseline_vasopressor"),
      " 和 ",
      get_pct(all_data, "baseline_crrt_rrt"),
      "。三队列在基线疾病严重程度、酸碱/氧合状态和治疗支持方面存在差异（Table 1；Table S1）。"
    )
  )

  text
}

result1_text <- manuscript_stats()

quote_powershell_path <- function(path) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "\\", mustWork = FALSE)), "'")
}

export_pptx_assets <- function(pptx_path, pdf_path, png_path) {
  if (.Platform$OS.type != "windows") {
    warning("Skipping Figure S1 PDF/PNG export: PowerPoint COM export is Windows-only.", call. = FALSE)
    return(invisible(FALSE))
  }

  pptx_path <- normalizePath(pptx_path, winslash = "\\", mustWork = TRUE)
  pdf_path <- normalizePath(pdf_path, winslash = "\\", mustWork = FALSE)
  png_path <- normalizePath(png_path, winslash = "\\", mustWork = FALSE)

  ps_script <- tempfile(fileext = ".ps1")
  writeLines(
    c(
      "$ErrorActionPreference = 'Stop'",
      paste0("$pptx = ", quote_powershell_path(pptx_path)),
      paste0("$pdf = ", quote_powershell_path(pdf_path)),
      paste0("$png = ", quote_powershell_path(png_path)),
      "$powerpoint = New-Object -ComObject PowerPoint.Application",
      "$presentation = $powerpoint.Presentations.Open($pptx, $true, $false, $false)",
      "try {",
      "  $presentation.SaveAs($pdf, 32)",
      "  $presentation.Slides.Item(1).Export($png, 'PNG', 3508, 2480) | Out-Null",
      "}",
      "finally {",
      "  $presentation.Close()",
      "  $powerpoint.Quit()",
      "  [System.Runtime.InteropServices.Marshal]::ReleaseComObject($presentation) | Out-Null",
      "  [System.Runtime.InteropServices.Marshal]::ReleaseComObject($powerpoint) | Out-Null",
      "  [GC]::Collect()",
      "  [GC]::WaitForPendingFinalizers()",
      "}"
    ),
    con = ps_script,
    useBytes = TRUE
  )

  status <- system2(
    "powershell",
    c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", shQuote(ps_script))
  )
  unlink(ps_script)
  if (!identical(status, 0L)) {
    warning("Failed to export Figure S1 PDF/PNG from PPTX.", call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# 导出结果 -------------------------------------------------------------------

write_csv(
  cohort_flow %>%
    select(
      "analysis_cohort",
      "cohort_label",
      "step_order",
      "step_id",
      "label",
      "n_before_step",
      "n_excluded_this_step",
      "n_remaining_after_step"
    ),
  file.path(TABLE_SUPP_DIR, "Figure_S1_cohort_selection_flow_data.csv"),
  na = ""
)

message("Release candidate: skipping nonessential PPTX/PDF/PNG flowchart export; retaining flowchart CSV data only.")

write_csv(table1, file.path(TABLE_MAIN_DIR, "Table_1_baseline_characteristics.csv"), na = "")
write.xlsx(table1, file.path(TABLE_MAIN_DIR, "Table_1_baseline_characteristics.xlsx"), overwrite = TRUE)

write_csv(table_s1, file.path(TABLE_SUPP_DIR, "Table_S1_missingness.csv"), na = "")
write.xlsx(table_s1, file.path(TABLE_SUPP_DIR, "Table_S1_missingness.xlsx"), overwrite = TRUE)

writeLines(
  result1_text,
  con = file.path(MANUSCRIPT_DIR, "simulated_results_result1_zh.md"),
  useBytes = TRUE
)

result1_qc <- tibble(
  item = c(
    "included_patients",
    "death_28d",
    "death_28d_pct",
    "observable_day_windows",
    "observable_day_windows_missing_twa_paco2"
  ),
  value = c(
    nrow(included_baseline),
    sum(included_baseline$death_28d %in% TRUE),
    round(sum(included_baseline$death_28d %in% TRUE) / nrow(included_baseline) * 100, 2),
    nrow(included_day),
    sum(is.na(included_day$twa_paco2))
  )
)

write_csv(result1_qc, file.path(TABLE_SUPP_DIR, "Result_1_qc_summary.csv"), na = "")

message("Done Result 1.")
message("Figure S1 flow data: ", file.path(TABLE_SUPP_DIR, "Figure_S1_cohort_selection_flow_data.csv"))
message("Table 1: ", file.path(TABLE_MAIN_DIR, "Table_1_baseline_characteristics.xlsx"))
message("Table S1: ", file.path(TABLE_SUPP_DIR, "Table_S1_missingness.xlsx"))
message("Result text: ", file.path(MANUSCRIPT_DIR, "simulated_results_result1_zh.md"))
