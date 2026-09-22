#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 05_build_report.R CONFIG_YAML RUN_DIR")
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = TRUE)
script_dir <- dirname(config_path)
source(file.path(script_dir, "R", "common.R"))

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)
activate_task_library(run_dir, cfg)
if (!identical(Sys.info()[["nodename"]], cfg$required_hostname)) {
  stop("Required host is ", cfg$required_hostname, "; observed ", Sys.info()[["nodename"]])
}
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("rmarkdown is unavailable in the selected SIF")

inventory_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
inventory_files <- inventory_files[file.info(inventory_files)$isdir %in% FALSE]
inventory <- data.frame(
  relative_path = substring(inventory_files, nchar(run_dir) + 2L),
  size_bytes = file.info(inventory_files)$size,
  modified_at = format(file.info(inventory_files)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_tsv(inventory, file.path(run_dir, "provenance", "run_manifest.tsv"))

old_env <- Sys.getenv(c("CARTKFLS_RUN_DIR", "CARTKFLS_SCRIPT_DIR"), unset = NA_character_)
on.exit({
  for (nm in names(old_env)) {
    if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
  }
}, add = TRUE)
Sys.setenv(CARTKFLS_RUN_DIR = run_dir, CARTKFLS_SCRIPT_DIR = script_dir)
output_path <- rmarkdown::render(
  input = file.path(script_dir, "report_template.Rmd"),
  output_file = "harmonization_technical_report.html",
  output_dir = file.path(run_dir, "report"),
  intermediates_dir = file.path(run_dir, "report", "intermediates"),
  clean = TRUE,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
if (!file.exists(output_path) || file.info(output_path)$size < 10000) stop("Report render did not produce a valid HTML file")
html_lines <- readLines(output_path, warn = FALSE)
embedded_png_count <- sum(lengths(regmatches(
  html_lines, gregexpr("data:image/png;base64", html_lines, fixed = TRUE)
)))
if (embedded_png_count < 4L) {
  stop("Report failed image embedding validation: observed ", embedded_png_count, " embedded PNGs; expected at least 4")
}
writeLines(c(
  paste0("completed_at=", timestamp()),
  paste0("hostname=", Sys.info()[["nodename"]]),
  paste0("report=", normalizePath(output_path)),
  paste0("report_size_bytes=", file.info(output_path)$size),
  paste0("embedded_png_count=", embedded_png_count)
), file.path(run_dir, "provenance", "report_complete.txt"))
log_message("Technical report complete: ", output_path)
