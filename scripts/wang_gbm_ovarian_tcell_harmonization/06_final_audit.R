#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 06_final_audit.R CONFIG_YAML RUN_DIR")
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = TRUE)
script_dir <- dirname(config_path)
source(file.path(script_dir, "R", "common.R"))

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)
activate_task_library(run_dir, cfg)
suppressPackageStartupMessages(library(data.table))
if (!identical(Sys.info()[["nodename"]], cfg$required_hostname)) {
  stop("Required host is ", cfg$required_hostname, "; observed ", Sys.info()[["nodename"]])
}

read_result <- function(...) fread(file.path(run_dir, ...), na.strings = "NA")
qc <- read_result("qc_tables", "harmonization_validation.tsv")
coverage <- rbindlist(list(
  read_result("qc_tables", "GBM_signature_gene_coverage.tsv"),
  read_result("qc_tables", "Ovarian_signature_gene_coverage.tsv")
))
counts <- read_result("qc_tables", "harmonized_annotation_counts.tsv")
patients <- read_result("sample_summaries", "harmonized_patient_summary.tsv")
best <- read_result("association_results", "best_parameter_row_audit.tsv")
effects <- read_result("association_results", "patient_chromosome_fitness_effects.tsv")
associations <- read_result("association_results", "chromosome_feature_spearman.tsv")
primary <- read_result("association_results", "primary_total_t_associations.tsv")
adjustments <- read_result("association_results", "one_feature_adjustment_models.tsv")
cross <- read_result("association_results", "cross_cancer_direction_comparison.tsv")

checks <- list()
add_check <- function(name, observed, expected, passed) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, observed = as.character(observed), expected = as.character(expected), passed = isTRUE(passed)
  )
}

sentinels <- c(
  "audit/audit_complete.txt", "provenance/harmonization_complete.txt",
  "provenance/qc_complete.txt", "provenance/fitness_association_complete.txt",
  "provenance/report_complete.txt", "provenance/pipeline_complete.txt",
  "provenance/postanalysis_refresh_complete.txt",
  "provenance/annotated_object_readback_complete.txt"
)
sentinel_present <- file.exists(file.path(run_dir, sentinels))
add_check("all_completion_sentinels_present", sum(sentinel_present), length(sentinels), all(sentinel_present))
add_check("all_harmonization_qc_checks_pass", sum(qc$passed %in% TRUE), nrow(qc), nrow(qc) == 6L && all(qc$passed %in% TRUE))
add_check("minimum_signature_gene_coverage", min(coverage$coverage_fraction), 1, min(coverage$coverage_fraction) == 1)

candidate_totals <- counts[, .(n_candidates = sum(N)), by = dataset]
add_check("GBM_candidate_cell_count", candidate_totals[dataset == "GBM", n_candidates], 924, candidate_totals[dataset == "GBM", n_candidates] == 924)
add_check("Ovarian_candidate_cell_count", candidate_totals[dataset == "Ovarian", n_candidates], 34527, candidate_totals[dataset == "Ovarian", n_candidates] == 34527)

patient_counts <- patients[, .(n_patients = uniqueN(patient_id)), by = dataset]
add_check("GBM_patient_summary_count", patient_counts[dataset == "GBM", n_patients], 20, patient_counts[dataset == "GBM", n_patients] == 20)
add_check("Ovarian_patient_summary_count", patient_counts[dataset == "Ovarian", n_patients], 11, patient_counts[dataset == "Ovarian", n_patients] == 11)

selected <- best[primary_patient_row %in% TRUE, .(n_selected_patients = uniqueN(PatientID)), by = dataset]
effect_expected <- sum(selected$n_selected_patients) * as.integer(cfg$fitness$chromosomes)
add_check("patient_chromosome_effect_row_count", nrow(effects), effect_expected, nrow(effects) == effect_expected)
add_check("feature_association_row_count", nrow(associations), 1452, nrow(associations) == 1452L)
add_check("primary_total_T_association_row_count", nrow(primary), 132, nrow(primary) == 132L)
add_check("adjustment_FDR_columns_present", paste(intersect(c("added_feature_fdr_bh", "adjusted_total_t_fdr_bh"), colnames(adjustments)), collapse = ";"), "added_feature_fdr_bh;adjusted_total_t_fdr_bh", all(c("added_feature_fdr_bh", "adjusted_total_t_fdr_bh") %in% colnames(adjustments)))

object_paths <- file.path(run_dir, "annotated_objects", c(
  "GBM_harmonized_tcell_annotations.Rds", "Ovarian_harmonized_tcell_annotations.Rds"
))
object_sizes <- file.info(object_paths)$size
add_check("two_nontrivial_harmonized_objects", paste(object_sizes, collapse = ";"), "both_gt_1GB", length(object_sizes) == 2L && all(object_sizes > 1e9))

report_path <- file.path(run_dir, "report", "harmonization_technical_report.html")
html_lines <- readLines(report_path, warn = FALSE)
embedded_png_count <- sum(lengths(regmatches(
  html_lines, gregexpr("data:image/png;base64", html_lines, fixed = TRUE)
)))
add_check("report_embedded_PNG_count", embedded_png_count, ">=4", embedded_png_count >= 4L)

audit <- rbindlist(checks)
fwrite(audit, file.path(run_dir, "provenance", "final_audit_summary.tsv"), sep = "\t", na = "NA")
if (any(!audit$passed)) stop("Final audit failed: ", paste(audit[passed == FALSE, check], collapse = ", "))

subtype_totals <- counts[, .(n_cells = sum(N)), by = .(dataset, harmonized_tcell_subtype)]
primary_fdr_hits <- primary[is.finite(fdr_bh) & fdr_bh < 0.05, .N]
feature_fdr_hits <- associations[is.finite(fdr_bh) & fdr_bh < 0.05, .N]
adjustment_fdr_hits <- adjustments[is.finite(added_feature_fdr_bh) & added_feature_fdr_bh < 0.05, .N]
opposite_programs <- cross[opposite_direction %in% TRUE & grepl("^mean_program_", feature)]
setorder(opposite_programs, -absolute_direction_difference)
strongest_opposite <- if (nrow(opposite_programs)) opposite_programs[1L] else data.table()

findings <- rbindlist(list(
  subtype_totals[, .(metric = paste(dataset, "subtype", harmonized_tcell_subtype, sep = ":"), value = as.character(n_cells), note = "cell count")],
  data.table(
    metric = c(
      "primary_total_T_FDR_hits", "all_feature_FDR_hits", "adjustment_added_feature_FDR_hits",
      "opposite_direction_program_rows", "strongest_opposite_program"
    ),
    value = c(
      primary_fdr_hits, feature_fdr_hits, adjustment_fdr_hits,
      nrow(opposite_programs),
      if (nrow(strongest_opposite)) paste0(
        "chr", strongest_opposite$chromosome, ":", strongest_opposite$chromosome_effect_metric,
        ":", strongest_opposite$feature, ":GBM=", signif(strongest_opposite$GBM, 4),
        ":Ovarian=", signif(strongest_opposite$Ovarian, 4)
      ) else "none"
    ),
    note = c(
      "BH FDR < 0.05", "BH FDR < 0.05; includes low-n perfect-rank correlations",
      "BH FDR < 0.05", "same feature and chromosome, opposite Spearman signs",
      "largest absolute cross-cancer direction difference"
    )
  )
), use.names = TRUE)
fwrite(findings, file.path(run_dir, "association_results", "final_key_findings.tsv"), sep = "\t", na = "NA")

inventory_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
inventory_files <- inventory_files[file.info(inventory_files)$isdir %in% FALSE]
inventory <- data.frame(
  relative_path = substring(inventory_files, nchar(run_dir) + 2L),
  size_bytes = file.info(inventory_files)$size,
  modified_at = format(file.info(inventory_files)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_tsv(inventory, file.path(run_dir, "provenance", "run_manifest.tsv"))
writeLines(c(
  paste0("completed_at=", timestamp()),
  paste0("hostname=", Sys.info()[["nodename"]]),
  paste0("checks_passed=", sum(audit$passed)),
  paste0("checks_total=", nrow(audit)),
  "final_audit_status=PASS"
), file.path(run_dir, "provenance", "final_audit_complete.txt"))
log_message("Final audit complete: PASS")
