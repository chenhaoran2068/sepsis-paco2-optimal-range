# Release-candidate runtime configuration.
# This file deliberately rejects implicit internal-project paths.

release_env <- function(name, default = NULL, required = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    if (required) {
      stop("Missing required environment variable: ", name, call. = FALSE)
    }
    return(default)
  }
  value
}

normalize_existing_dir <- function(path, label) {
  if (!dir.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

is_within <- function(path, parent) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent <- normalizePath(parent, winslash = "/", mustWork = FALSE)
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

release_paths <- function() {
  release_root <- normalize_existing_dir(
    release_env("SEPSIS_PACO2_RELEASE_ROOT", required = TRUE),
    "SEPSIS_PACO2_RELEASE_ROOT"
  )
  requested_mode <- release_env("SEPSIS_PACO2_RELEASE_MODE", default = "demo")
  if (!identical(requested_mode, "demo")) {
    stop("This public reproducibility package supports simulated demo data only.", call. = FALSE)
  }
  mode <- "demo"

  run_id <- release_env("SEPSIS_PACO2_RUN_ID", default = "demo_validation_v0_1_0")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", run_id)) {
    stop("SEPSIS_PACO2_RUN_ID contains unsupported characters.", call. = FALSE)
  }

  default_input <- file.path(release_root, "demo", "input")
  if (nzchar(Sys.getenv("SEPSIS_PACO2_ANALYSIS_DATA_ROOT", unset = ""))) {
    stop("This public package does not permit an external analysis-data root.", call. = FALSE)
  }
  input_root <- normalize_existing_dir(default_input, "simulated demo input")

  if (nzchar(Sys.getenv("SEPSIS_PACO2_RESULT_ROOT", unset = ""))) {
    stop("This public package does not permit an external result root.", call. = FALSE)
  }
  default_result <- file.path(release_root, "output", "demo_runs", run_id)
  result_root <- normalizePath(default_result, winslash = "/", mustWork = FALSE)
  expected_parent <- file.path(release_root, "output", "demo_runs")
  if (!is_within(result_root, expected_parent)) {
    stop("Result root must be below ", expected_parent, ".", call. = FALSE)
  }
  if (grepl("locked_baseline|/reruns/", result_root, ignore.case = TRUE)) {
    stop("The demonstration must not write to internal analysis or release-control paths.", call. = FALSE)
  }

  expected_script_root <- file.path(release_root, "code", "R")
  requested_script_root <- release_env("SEPSIS_PACO2_SCRIPT_ROOT", default = expected_script_root)
  script_root <- normalize_existing_dir(requested_script_root, "public script root")
  if (!identical(script_root, normalizePath(expected_script_root, winslash = "/", mustWork = TRUE))) {
    stop("This public package does not permit a replacement script root.", call. = FALSE)
  }

  seed <- suppressWarnings(as.integer(release_env("SEPSIS_PACO2_RANDOM_SEED", default = "20260908")))
  nthreads <- suppressWarnings(as.integer(release_env("SEPSIS_PACO2_NTHREADS", default = "1")))
  if (is.na(seed) || is.na(nthreads) || nthreads < 1L) {
    stop("Random seed and thread count must be valid positive integers.", call. = FALSE)
  }

  list(
    release_root = release_root,
    mode = mode,
    run_id = run_id,
    input_root = input_root,
    result_root = result_root,
    pamm_root = file.path(result_root, "intermediate", "pamm"),
    script_root = script_root,
    seed = seed,
    nthreads = nthreads
  )
}
