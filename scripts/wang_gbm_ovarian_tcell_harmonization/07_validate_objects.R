#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 07_validate_objects.R CONFIG_YAML RUN_DIR")
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = TRUE)
script_dir <- dirname(config_path)
source(file.path(script_dir, "R", "common.R"))

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)
activate_task_library(run_dir, cfg)
suppressPackageStartupMessages(library(Seurat))
if (!identical(Sys.info()[["nodename"]], cfg$required_hostname)) {
  stop("Required host is ", cfg$required_hostname, "; observed ", Sys.info()[["nodename"]])
}

programs <- c(
  "cytotoxicity", "activation", "naive_memory", "exhaustion_dysfunction",
  "interferon_response", "stress"
)
required_columns <- c(
  "harmonized_tnk_candidate", "harmonized_lineage", "harmonized_tcell_subtype",
  "harmonized_subtype_confidence", "harmonized_uncertainty_reason",
  paste0("harmonized_score_", programs)
)
contracts <- data.frame(
  dataset = c("GBM", "Ovarian"),
  filename = c("GBM_harmonized_tcell_annotations.Rds", "Ovarian_harmonized_tcell_annotations.Rds"),
  expected_cells = c(106700L, 91654L),
  expected_candidates = c(924L, 34527L),
  stringsAsFactors = FALSE
)

results <- lapply(seq_len(nrow(contracts)), function(i) {
  dataset <- contracts$dataset[[i]]
  path <- file.path(run_dir, "annotated_objects", contracts$filename[[i]])
  log_message("Readback validating ", dataset, " object: ", path)
  obj <- readRDS(path)
  meta <- obj@meta.data
  missing_columns <- setdiff(required_columns, colnames(meta))
  candidates <- meta$harmonized_tnk_candidate %in% TRUE
  complete_programs <- if (!length(missing_columns)) {
    all(stats::complete.cases(meta[candidates, paste0("harmonized_score_", programs), drop = FALSE]))
  } else {
    FALSE
  }
  row <- data.frame(
    dataset = dataset,
    object_path = path,
    object_size_bytes = file.info(path)$size,
    object_class = paste(class(obj), collapse = ";"),
    observed_cells = ncol(obj),
    expected_cells = contracts$expected_cells[[i]],
    observed_candidates = sum(candidates),
    expected_candidates = contracts$expected_candidates[[i]],
    missing_required_columns = paste(missing_columns, collapse = ";"),
    all_candidate_program_scores_complete = complete_programs,
    passed = inherits(obj, "Seurat") && ncol(obj) == contracts$expected_cells[[i]] &&
      sum(candidates) == contracts$expected_candidates[[i]] && !length(missing_columns) && complete_programs,
    stringsAsFactors = FALSE
  )
  rm(obj, meta)
  invisible(gc(verbose = FALSE))
  row
})
results <- do.call(rbind, results)
write_tsv(results, file.path(run_dir, "provenance", "annotated_object_readback_validation.tsv"))
if (any(!results$passed)) stop("Annotated object readback validation failed")
writeLines(c(
  paste0("completed_at=", timestamp()),
  paste0("hostname=", Sys.info()[["nodename"]]),
  paste0("objects_validated=", nrow(results)),
  "annotated_object_readback_status=PASS"
), file.path(run_dir, "provenance", "annotated_object_readback_complete.txt"))
log_message("Annotated object readback validation complete: PASS")
