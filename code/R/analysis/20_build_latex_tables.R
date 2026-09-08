options(stringsAsFactors = FALSE)

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
TABLE_MAIN_DIR <- file.path(RESULT_ROOT, "tables", "main")
TABLE_SUPP_DIR <- file.path(RESULT_ROOT, "tables", "supplementary")
LATEX_TABLE_DIR <- file.path(RESULT_ROOT, "tables", "latex")
dir.create(LATEX_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub(">=", "$\\ge$", x, fixed = TRUE)
  x <- gsub("<=", "$\\le$", x, fixed = TRUE)
  x <- gsub("<", "$<$", x, fixed = TRUE)
  x <- gsub(">", "$>$", x, fixed = TRUE)
  x <- gsub("PaCO2", "PaCO$_2$", x, fixed = TRUE)
  x <- gsub("TWA-PaCO$_2$", "TWA-PaCO$_2$", x, fixed = TRUE)
  x <- gsub("HCO3", "HCO$_3^-$", x, fixed = TRUE)
  x <- gsub("PaO2/FiO2", "PaO$_2$/FiO$_2$", x, fixed = TRUE)
  x <- gsub("Day 1-7", "Days 1--7", x, fixed = TRUE)
  x <- gsub("1-7", "1--7", x, fixed = TRUE)
  x
}

write_tabular <- function(data, path, col_spec, header = names(data), escape_header = TRUE) {
  escaped <- as.data.frame(lapply(data, latex_escape), check.names = FALSE)
  header_line <- if (escape_header) {
    paste(latex_escape(header), collapse = " & ")
  } else {
    paste(header, collapse = " & ")
  }
  lines <- c(
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    header_line,
    "\\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(escaped))) {
    lines <- c(lines, paste(unlist(escaped[i, ], use.names = FALSE), collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, path, useBytes = TRUE)
}

write_longtable <- function(data, path, col_spec, header = names(data), escape_header = TRUE) {
  escaped <- as.data.frame(lapply(data, latex_escape), check.names = FALSE)
  header_line <- if (escape_header) {
    paste(latex_escape(header), collapse = " & ")
  } else {
    paste(header, collapse = " & ")
  }
  lines <- c(
    paste0("\\begin{longtable}{", col_spec, "}"),
    "\\toprule",
    header_line,
    "\\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    header_line,
    "\\\\",
    "\\midrule",
    "\\endhead"
  )
  for (i in seq_len(nrow(escaped))) {
    lines <- c(lines, paste(unlist(escaped[i, ], use.names = FALSE), collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{longtable}")
  writeLines(lines, path, useBytes = TRUE)
}

read_csv <- function(path) read.csv(path, check.names = FALSE)

table1 <- read_csv(file.path(TABLE_MAIN_DIR, "Table_1_baseline_characteristics.csv"))
table1 <- table1[, c("Section", "Characteristic", "Display", "Overall", "MIMIC-IV", "AmsterdamUMCdb", "Chinese")]
names(table1)[3] <- "Statistic"
write_tabular(
  table1,
  file.path(LATEX_TABLE_DIR, "table1.tex"),
  "@{}p{0.16\\linewidth}p{0.24\\linewidth}p{0.11\\linewidth}p{0.13\\linewidth}p{0.13\\linewidth}p{0.13\\linewidth}p{0.13\\linewidth}@{}"
)

table2 <- read_csv(file.path(TABLE_MAIN_DIR, "Table_2_primary_landmark_categories.csv"))
write_tabular(table2, file.path(LATEX_TABLE_DIR, "table2.tex"), "@{}llllllll@{}")

table3 <- read_csv(file.path(TABLE_MAIN_DIR, "Table_3_high_risk_paco2_burden.csv"))
write_tabular(
  table3,
  file.path(LATEX_TABLE_DIR, "table3.tex"),
  "@{}p{0.14\\linewidth}p{0.09\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.08\\linewidth}p{0.08\\linewidth}p{0.08\\linewidth}p{0.06\\linewidth}p{0.06\\linewidth}p{0.08\\linewidth}p{0.08\\linewidth}p{0.08\\linewidth}@{}",
  header = c(
    "Exposure metric",
    "\\shortstack{Category\\\\or contrast}",
    "Reference",
    "\\shortstack{Landmark\\\\days}",
    "\\shortstack{N at\\\\landmark}",
    "\\shortstack{N in\\\\category}",
    "\\shortstack{Events in\\\\category}",
    "\\shortstack{N days\\\\estimated}",
    "Adjusted HR",
    "95\\% CI",
    "P value",
    "\\shortstack{P for\\\\trend}"
  ),
  escape_header = FALSE
)

supp_specs <- list(
  table_s1 = list(file = "Table_S1_missingness.csv", spec = "@{}lllll@{}", type = "tabular"),
  table_s2 = list(file = "Table_S2_pamm_adjustment_hierarchy.csv", spec = "@{}lllllllllllll@{}", type = "tabular",
                  header = c("Model", "Set", "Patients", "PED rows", "Events", "Nadir", "5\\% range", "Smooth P",
                             "\\shortstack{HR 35\\\\vs nadir}", "\\shortstack{HR 50\\\\vs nadir}",
                             "\\shortstack{HR 60\\\\vs nadir}", "\\shortstack{HR 50\\\\vs 40}", "AIC")),
  table_s3 = list(file = "Table_S3_fine_paco2_category_sensitivity_landmark.csv", spec = "@{}llllllll@{}", type = "tabular"),
  table_s4 = list(file = "Table_S4_paco2_category_distribution_event_rates.csv", spec = "@{}llllllll@{}", type = "longtable"),
  table_s5 = list(file = "Table_S5_secondary_exposure_definitions.csv", spec = "@{}>{\\raggedright\\arraybackslash}p{1.30in}>{\\raggedright\\arraybackslash}p{1.20in}>{\\raggedright\\arraybackslash}p{0.95in}>{\\raggedright\\arraybackslash}p{4.25in}@{}", type = "tabular"),
  table_s6 = list(file = "Table_S6_patient_level_high_risk_burden_support.csv", spec = "@{}>{\\raggedright\\arraybackslash}p{0.90in}>{\\raggedright\\arraybackslash}p{1.35in}>{\\raggedright\\arraybackslash}p{1.20in}rr>{\\raggedright\\arraybackslash}p{0.75in}>{\\raggedright\\arraybackslash}p{1.55in}>{\\raggedright\\arraybackslash}p{1.20in}@{}", type = "tabular"),
  table_s7 = list(file = "Table_S7_subgroup_estimates.csv", spec = "@{}llllllll@{}", type = "tabular"),
  table_s8 = list(file = "Table_S8_ph_adjusted_sensitivity.csv", spec = "@{}llllllllll@{}", type = "tabular",
                  header = c("Model", "Patients", "PED rows", "Events", "Nadir", "5\\% range",
                             "\\shortstack{HR at\\\\35}", "\\shortstack{HR at\\\\50}", "\\shortstack{HR at\\\\60}", "AIC"))
)

for (nm in names(supp_specs)) {
  spec <- supp_specs[[nm]]
  dat <- read_csv(file.path(TABLE_SUPP_DIR, spec$file))
  out <- file.path(LATEX_TABLE_DIR, paste0(nm, ".tex"))
  if (identical(spec$type, "longtable")) {
    write_longtable(dat, out, spec$spec)
  } else {
    if (!is.null(spec$header)) {
      write_tabular(dat, out, spec$spec, header = spec$header, escape_header = FALSE)
    } else {
      write_tabular(dat, out, spec$spec)
    }
  }
}

message("LaTeX table files written to: ", LATEX_TABLE_DIR)
