#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: test_postharmonization_contract.R RUN_DIR")
run_dir <- normalizePath(args[[1L]], mustWork = TRUE)
suppressPackageStartupMessages(library(data.table))

programs <- c(
  "cytotoxicity", "activation", "naive_memory", "exhaustion_dysfunction",
  "interferon_response", "stress"
)
required_common <- c(
  "harmonized_tnk_candidate", "harmonized_lineage", "harmonized_tcell_subtype",
  "harmonized_subtype_confidence", "harmonized_uncertainty_reason",
  paste0("harmonized_score_", programs)
)

gbm_path <- file.path(run_dir, "cell_metadata", "GBM_harmonized_cell_metadata.tsv.gz")
ovarian_path <- file.path(run_dir, "cell_metadata", "Ovarian_harmonized_cell_metadata.tsv.gz")
stopifnot(file.exists(gbm_path), file.exists(ovarian_path))
gbm <- fread(gbm_path, nrows = 100)
ovarian <- fread(ovarian_path, nrows = 100)
stopifnot(
  !length(setdiff(c(required_common, "PatientID", "Sample"), colnames(gbm))),
  !length(setdiff(c(required_common, "patient_id", "sample_id", "cell_type_subtype_final"), colnames(ovarian)))
)
state_column <- intersect(c("cell_state_final", "cell_state_program_final", "cell_state_stage1"), colnames(ovarian))
stopifnot(length(state_column) >= 1L)
state_column <- state_column[[1L]]
crosswalk <- ovarian[, .N, by = c(
  "cell_type_subtype_final", state_column,
  "harmonized_tcell_subtype", "harmonized_state_primary"
)]
setnames(crosswalk, c("cell_type_subtype_final", state_column), c("original_subtype", "original_state"))
stopifnot(all(c("original_subtype", "original_state", "N") %in% colnames(crosswalk)))
cat("postharmonization_contract=PASS\n")
