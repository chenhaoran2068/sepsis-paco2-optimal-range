# Run the complete simulated-data demonstration without touching study data.

args <- commandArgs(trailingOnly = TRUE)
run_id <- "demo_run_001"
if (length(args) > 0L) {
  if (length(args) != 2L || args[[1]] != "--run-id" || !nzchar(args[[2]])) {
    stop("Usage: Rscript code/run_demo.R --run-id <new-run-id>", call. = FALSE)
  }
  run_id <- args[[2]]
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Unable to determine the run_demo.R location.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
release_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
script_root <- file.path(release_root, "code", "R")
source(file.path(script_root, "00_functions", "release_config.R"))

Sys.setenv(
  SEPSIS_PACO2_RELEASE_ROOT = release_root,
  SEPSIS_PACO2_RELEASE_MODE = "demo",
  SEPSIS_PACO2_RUN_ID = run_id,
  SEPSIS_PACO2_RANDOM_SEED = "20260908",
  SEPSIS_PACO2_NTHREADS = "1",
  SEPSIS_PACO2_SCRIPT_ROOT = script_root
)
dir.create(file.path(release_root, "demo", "input"), recursive = TRUE, showWarnings = FALSE)
paths <- release_paths()

if (dir.exists(paths$result_root)) {
  stop("Refusing to overwrite existing demo run: ", paths$result_root, call. = FALSE)
}
longest_expected_output <- c(
  "figures/supplementary/Figure_S3_pamm_adjustment_hierarchy_model_qc.png",
  "figures/supplementary/Figure_S5_cohort_day_specific_landmark_summary.png"
)
if (.Platform$OS.type == "windows" && max(nchar(file.path(paths$result_root, longest_expected_output))) >= 259L) {
  stop(
    "The configured demo output path is too long for reliable Windows figure output. ",
    "Use a shorter repository location or a shorter --run-id.",
    call. = FALSE
  )
}
dir.create(paths$result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$result_root, "logs"), recursive = TRUE, showWarnings = FALSE)

demo_label <- c(
  "# Simulated Data Demonstration",
  "",
  "This directory contains synthetic-data pipeline output only.",
  "It is not derived from study participants and must not be cited as manuscript results."
)
writeLines(demo_label, file.path(paths$result_root, "00_SIMULATED_DEMO_ONLY.md"))

rscript_name <- if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript"
rscript <- file.path(R.home("bin"), rscript_name)
if (!file.exists(rscript)) {
  stop("Unable to locate Rscript for child analysis stages.", call. = FALSE)
}

run_stage <- function(name, script) {
  log <- file.path(paths$result_root, "logs", paste0(name, ".log"))
  status <- system2(rscript, args = shQuote(script), stdout = log, stderr = log)
  if (!identical(status, 0L)) {
    stop("Stage failed: ", name, ". See ", log, call. = FALSE)
  }
  message("Completed: ", name)
}

run_stage("00_generate_demo_data", file.path(script_root, "generate_demo_data.R"))
run_stage("01_prepare_pamm_ped_data", file.path(script_root, "analysis", "01_prepare_pamm_ped_data.R"))
run_stage("10_result1_population_baseline", file.path(script_root, "analysis", "10_result1_population_baseline.R"))
run_stage("11_result2_primary_pamm", file.path(script_root, "analysis", "11_result2_primary_pamm.R"))
run_stage("13_result3_landmark_burden", file.path(script_root, "analysis", "13_result3_landmark_burden.R"))
run_stage("14_result4_subgroups_ph_sensitivity", file.path(script_root, "analysis", "14_result4_subgroups_ph_sensitivity.R"))
run_stage("20_build_latex_tables", file.path(script_root, "analysis", "20_build_latex_tables.R"))

for (directory in c("figures", "tables", "manuscript", "intermediate")) {
  label_path <- file.path(paths$result_root, directory, "00_SIMULATED_DEMO_ONLY.md")
  dir.create(dirname(label_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(demo_label, label_path)
}

run_stage("11_validate_demo_run", file.path(release_root, "tests", "validate_demo_run.R"))
message("Simulated-data demo completed: ", paths$result_root)
