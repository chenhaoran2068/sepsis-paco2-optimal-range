# Validate the public simulated-data demonstration without inspecting study data.

release_root <- Sys.getenv("SEPSIS_PACO2_RELEASE_ROOT", unset = "")
script_root <- Sys.getenv("SEPSIS_PACO2_SCRIPT_ROOT", unset = "")
if (!nzchar(release_root) || !nzchar(script_root)) {
  stop("Set release configuration through code/run_demo.R.", call. = FALSE)
}
source(file.path(script_root, "00_functions", "release_config.R"))
paths <- release_paths()

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
assert_file <- function(path) {
  assert_true(file.exists(path), paste0("Missing required file: ", path))
  assert_true(file.info(path)$size > 0L, paste0("Empty required file: ", path))
}

assert_true(paths$mode == "demo", "QA is restricted to demo mode.")
assert_true(is_within(paths$result_root, file.path(paths$release_root, "output", "demo_runs")), "Unexpected result root.")
assert_file(file.path(paths$result_root, "00_SIMULATED_DEMO_ONLY.md"))

if (!requireNamespace("arrow", quietly = TRUE)) {
  stop("Package 'arrow' is required for demo QA.", call. = FALSE)
}

input_root <- paths$input_root
baseline_path <- file.path(input_root, "pooled", "all_baseline_outcome.parquet")
day_path <- file.path(input_root, "pooled", "all_day_long.parquet")
raw_path <- file.path(input_root, "pooled", "all_raw_paco2.parquet")
for (path in c(baseline_path, day_path, raw_path)) assert_file(path)

baseline <- as.data.frame(arrow::read_parquet(baseline_path))
day_long <- as.data.frame(arrow::read_parquet(day_path))
raw_paco2 <- as.data.frame(arrow::read_parquet(raw_path))
required_cohorts <- c("MIMIC", "AmsterdamUMCdb", "Chinese cohort")
assert_true(nrow(baseline) == 1500L, "Synthetic baseline dataset must contain 1,500 stays.")
assert_true(all(table(baseline$analysis_cohort)[required_cohorts] == 500L), "Synthetic cohort counts are invalid.")
assert_true(nrow(day_long) == 10500L, "Synthetic longitudinal dataset must contain 10,500 rows.")
assert_true(nrow(raw_paco2) == 10500L, "Synthetic PaCO2 dataset must contain 10,500 rows.")
assert_true(all(day_long$icu_day %in% 1:7), "Synthetic ICU days must be in 1-7.")
assert_true(!anyNA(baseline$global_stay_id), "Synthetic stay identifiers must be non-missing.")

contract_path <- file.path(paths$release_root, "expected", "output-contract.csv")
assert_file(contract_path)
contract <- read.csv(contract_path, check.names = FALSE, stringsAsFactors = FALSE)
assert_true(all(c("relative_path", "role") %in% names(contract)), "Output contract is invalid.")
required_outputs <- contract$relative_path
assert_true(length(required_outputs) > 0L && !anyDuplicated(required_outputs), "Output contract must list unique outputs.")
for (relative_path in required_outputs) assert_file(file.path(paths$result_root, relative_path))

table1 <- read.csv(file.path(paths$result_root, "tables", "main", "Table_1_baseline_characteristics.csv"), check.names = FALSE)
table2 <- read.csv(file.path(paths$result_root, "tables", "main", "Table_2_primary_landmark_categories.csv"), check.names = FALSE)
table3 <- read.csv(file.path(paths$result_root, "tables", "main", "Table_3_high_risk_paco2_burden.csv"), check.names = FALSE)
assert_true(all(c("Characteristic", "Overall", "MIMIC-IV", "AmsterdamUMCdb", "Chinese") %in% names(table1)), "Table 1 structure is invalid.")
assert_true(all(c("Analysis set", "Exposure category", "Adjusted HR", "95% CI") %in% names(table2)), "Table 2 structure is invalid.")
assert_true(all(c("Exposure metric", "Category or contrast", "Adjusted HR", "95% CI") %in% names(table3)), "Table 3 structure is invalid.")

public_surface <- c(
  file.path(paths$release_root, "README.md"),
  list.files(file.path(paths$release_root, "code"), pattern = "\\.(R|md|yml)$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(paths$release_root, "demo"), pattern = "\\.(R|md|csv)$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(paths$release_root, "expected"), pattern = "\\.(R|md|csv)$", recursive = TRUE, full.names = TRUE)
)
for (path in public_surface[file.exists(public_surface)]) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  has_windows_path <- grepl("(?<![[:alnum:]])[A-Za-z]:[/\\\\]", lines, perl = TRUE)
  assert_true(!any(has_windows_path), paste0("Absolute Windows path in public surface: ", path))
}

all_outputs <- list.files(paths$result_root, recursive = TRUE, full.names = TRUE)
all_outputs <- all_outputs[file.info(all_outputs)$isdir == FALSE]
result_prefix <- paste0(normalizePath(paths$result_root, winslash = "/"), "/")
inventory <- data.frame(
  relative_path = sub(result_prefix, "", normalizePath(all_outputs, winslash = "/"), fixed = TRUE),
  bytes = file.info(all_outputs)$size,
  md5 = unname(tools::md5sum(all_outputs)),
  stringsAsFactors = FALSE
)
qa_root <- file.path(paths$result_root, "qa")
dir.create(qa_root, recursive = TRUE, showWarnings = FALSE)
write.csv(inventory, file.path(qa_root, "demo_output_inventory.csv"), row.names = FALSE, na = "")
summary <- data.frame(
  status = "PASS",
  run_id = paths$run_id,
  mode = paths$mode,
  seed = paths$seed,
  baseline_rows = nrow(baseline),
  day_long_rows = nrow(day_long),
  raw_paco2_rows = nrow(raw_paco2),
  required_output_count = length(required_outputs),
  generated_file_count = nrow(inventory),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(qa_root, "demo_qa_summary.csv"), row.names = FALSE, na = "")
message("PASS: simulated-data demo QA completed for ", paths$run_id)
