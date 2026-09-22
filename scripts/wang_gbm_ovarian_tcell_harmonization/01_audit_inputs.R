#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_audit_inputs.R CONFIG_YAML RUN_DIR")
}
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = FALSE)
source(file.path(dirname(config_path), "R", "common.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)

host <- Sys.info()[["nodename"]]
if (!identical(host, cfg$required_hostname)) {
  stop("Refusing to run on host ", host, "; required host is ", cfg$required_hostname)
}

all_objects <- list()
all_columns <- list()
all_values <- list()

for (dataset in names(cfg$inputs)) {
  source_path <- cfg$inputs[[dataset]]
  if (!file.exists(source_path)) stop(dataset, " input not found: ", source_path)
  log_message("Reading ", dataset, " object: ", source_path)
  obj <- readRDS(source_path)
  if (!inherits(obj, "Seurat")) stop(dataset, " input is not a Seurat object")

  all_objects[[dataset]] <- object_profile(obj, dataset, source_path)
  p <- metadata_profile(obj@meta.data)
  p$dataset <- dataset
  p <- p[, c("dataset", setdiff(colnames(p), "dataset")), drop = FALSE]
  all_columns[[dataset]] <- p
  all_values[[dataset]] <- metadata_value_counts(obj@meta.data, dataset)

  saveRDS(
    obj@meta.data,
    file.path(run_dir, "audit", paste0(dataset, "_metadata_only.rds")),
    compress = TRUE
  )
  rm(obj)
  invisible(gc(verbose = FALSE))
  log_message("Completed metadata audit for ", dataset)
}

write_tsv(do.call(rbind, all_objects), file.path(run_dir, "audit", "input_object_profile.tsv"))
write_tsv(do.call(rbind, all_columns), file.path(run_dir, "audit", "metadata_column_profile.tsv"))
write_tsv(do.call(rbind, all_values), file.path(run_dir, "audit", "metadata_value_counts.tsv"))

writeLines(
  c(
    paste0("audit_completed_at=", timestamp()),
    paste0("hostname=", Sys.info()[["nodename"]]),
    paste0("R_version=", R.version.string),
    paste0("Seurat_version=", as.character(utils::packageVersion("Seurat")))
  ),
  file.path(run_dir, "audit", "audit_complete.txt")
)
log_message("Input audit complete")
