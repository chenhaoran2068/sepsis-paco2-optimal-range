# Restore the R package environment recorded for this private release candidate.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Unable to determine install_dependencies.R location.", call. = FALSE)
}
code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
release_root <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = TRUE)
lockfile <- file.path(release_root, "renv.lock")

if (!file.exists(lockfile)) {
  stop("Missing renv.lock: ", lockfile, call. = FALSE)
}
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::restore(project = release_root, lockfile = lockfile, prompt = FALSE)
message("R package environment restored for: ", release_root)
