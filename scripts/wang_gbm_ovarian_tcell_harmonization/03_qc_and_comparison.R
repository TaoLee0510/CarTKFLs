#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 03_qc_and_comparison.R CONFIG_YAML RUN_DIR")
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = TRUE)
script_dir <- dirname(config_path)
source(file.path(script_dir, "R", "common.R"))

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)
activate_task_library(run_dir, cfg)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!identical(Sys.info()[["nodename"]], cfg$required_hostname)) {
  stop("Required host is ", cfg$required_hostname, "; observed ", Sys.info()[["nodename"]])
}

program_signatures <- c(
  "cytotoxicity", "activation", "naive_memory", "exhaustion_dysfunction",
  "interferon_response", "stress"
)
subtype_levels <- c("CD8", "Conventional_CD4", "Treg", "Unresolved_T")

read_dataset_metadata <- function(dataset) {
  path <- file.path(run_dir, "cell_metadata", paste0(dataset, "_harmonized_cell_metadata.tsv.gz"))
  if (!file.exists(path)) stop("Missing harmonized metadata: ", path)
  x <- data.table::fread(path, sep = "\t", na.strings = c("NA", ""))
  if (dataset == "GBM") {
    x[, `:=`(patient_id_harmonized = as.character(PatientID), sample_id_harmonized = as.character(Sample))]
  } else {
    x[, `:=`(patient_id_harmonized = as.character(patient_id), sample_id_harmonized = as.character(sample_id))]
  }
  x
}

summarize_groups <- function(dt, group_columns, summary_level) {
  group_columns <- c("dataset", group_columns)
  base <- dt[, .(
    n_all_cells = .N,
    n_t_cells = sum(harmonized_lineage == "T_cell", na.rm = TRUE),
    n_nk_cells = sum(harmonized_lineage == "NK", na.rm = TRUE),
    n_tnk_candidates = sum(harmonized_tnk_candidate %in% TRUE, na.rm = TRUE)
  ), by = group_columns]
  base[, `:=`(
    total_t_fraction_all_cells = n_t_cells / n_all_cells,
    nk_fraction_all_cells = n_nk_cells / n_all_cells,
    summary_level = summary_level
  )]

  t_only <- dt[harmonized_lineage == "T_cell"]
  subtype <- t_only[, .N, by = c(group_columns, "harmonized_tcell_subtype")]
  if (nrow(subtype)) {
    subtype_wide <- data.table::dcast(
      subtype, stats::as.formula(paste(paste(group_columns, collapse = " + "), "~ harmonized_tcell_subtype")),
      value.var = "N", fill = 0
    )
    present <- intersect(subtype_levels, colnames(subtype_wide))
    data.table::setnames(subtype_wide, present, paste0("n_subtype_", present))
    base <- merge(base, subtype_wide, by = group_columns, all.x = TRUE, sort = FALSE)
  }
  for (level in subtype_levels) {
    n_col <- paste0("n_subtype_", level)
    p_col <- paste0("prop_within_t_", level)
    if (!n_col %in% colnames(base)) base[, (n_col) := 0L]
    base[is.na(get(n_col)), (n_col) := 0L]
    base[, (p_col) := ifelse(n_t_cells > 0, get(n_col) / n_t_cells, NA_real_)]
  }

  score_columns <- paste0("harmonized_score_", program_signatures)
  missing_scores <- setdiff(score_columns, colnames(t_only))
  if (length(missing_scores)) stop("Missing program scores: ", paste(missing_scores, collapse = ", "))
  long <- data.table::melt(
    t_only, id.vars = group_columns, measure.vars = score_columns,
    variable.name = "program", value.name = "score", variable.factor = FALSE
  )
  long[, program := sub("^harmonized_score_", "", program)]
  program_stats <- long[, .(
    mean_score = mean(score, na.rm = TRUE),
    median_score = stats::median(score, na.rm = TRUE),
    n_scored = sum(is.finite(score))
  ), by = c(group_columns, "program")]
  mean_wide <- data.table::dcast(
    program_stats, stats::as.formula(paste(paste(group_columns, collapse = " + "), "~ program")),
    value.var = "mean_score"
  )
  median_wide <- data.table::dcast(
    program_stats, stats::as.formula(paste(paste(group_columns, collapse = " + "), "~ program")),
    value.var = "median_score"
  )
  program_cols <- setdiff(colnames(mean_wide), group_columns)
  data.table::setnames(mean_wide, program_cols, paste0("mean_program_", program_cols))
  program_cols <- setdiff(colnames(median_wide), group_columns)
  data.table::setnames(median_wide, program_cols, paste0("median_program_", program_cols))
  base <- merge(base, mean_wide, by = group_columns, all.x = TRUE, sort = FALSE)
  merge(base, median_wide, by = group_columns, all.x = TRUE, sort = FALSE)
}

log_message("Reading harmonized cell metadata")
gbm <- read_dataset_metadata("GBM")
ovarian <- read_dataset_metadata("Ovarian")
all_meta <- data.table::rbindlist(list(gbm, ovarian), use.names = TRUE, fill = TRUE)

patient_summary <- summarize_groups(all_meta, "patient_id_harmonized", "patient")
data.table::setnames(patient_summary, "patient_id_harmonized", "patient_id")
sample_summary <- summarize_groups(
  all_meta, c("patient_id_harmonized", "sample_id_harmonized"), "sample"
)
data.table::setnames(
  sample_summary, c("patient_id_harmonized", "sample_id_harmonized"), c("patient_id", "sample_id")
)
data.table::fwrite(patient_summary, file.path(run_dir, "sample_summaries", "harmonized_patient_summary.tsv"), sep = "\t", na = "NA")
data.table::fwrite(sample_summary, file.path(run_dir, "sample_summaries", "harmonized_sample_summary.tsv"), sep = "\t", na = "NA")

count_table <- all_meta[harmonized_tnk_candidate %in% TRUE, .N, by = .(
  dataset, harmonized_lineage, harmonized_tcell_subtype,
  harmonized_subtype_confidence, harmonized_state_primary
)]
data.table::setorder(count_table, dataset, harmonized_lineage, harmonized_tcell_subtype, -N)
data.table::fwrite(count_table, file.path(run_dir, "qc_tables", "harmonized_annotation_counts.tsv"), sep = "\t", na = "NA")

reason_table <- all_meta[harmonized_tnk_candidate %in% TRUE, .N, by = .(
  dataset, harmonized_uncertainty_reason
)]
data.table::fwrite(reason_table, file.path(run_dir, "qc_tables", "uncertainty_reason_counts.tsv"), sep = "\t", na = "NA")

ovarian_state_column <- intersect(
  c("cell_state_final", "cell_state_program_final", "cell_state_stage1"),
  colnames(ovarian)
)
if (!length(ovarian_state_column)) stop("Ovarian metadata lacks an original state annotation column")
ovarian_state_column <- ovarian_state_column[[1L]]
ovarian_crosswalk <- ovarian[harmonized_tnk_candidate %in% TRUE, .N, by = c(
  "cell_type_subtype_final", ovarian_state_column,
  "harmonized_tcell_subtype", "harmonized_state_primary"
)]
data.table::setnames(
  ovarian_crosswalk,
  c("cell_type_subtype_final", ovarian_state_column),
  c("original_subtype", "original_state")
)
data.table::fwrite(ovarian_crosswalk, file.path(run_dir, "qc_tables", "Ovarian_original_to_harmonized_crosswalk.tsv"), sep = "\t", na = "NA")

qc <- data.table::rbindlist(list(
  all_meta[, .(
    check = "candidate_cells_have_subtype",
    passed = all(!harmonized_tnk_candidate | !is.na(harmonized_tcell_subtype)),
    observed = sum(harmonized_tnk_candidate & is.na(harmonized_tcell_subtype)), expected = 0L
  ), by = dataset],
  all_meta[, .(
    check = "candidate_cells_have_all_six_program_scores",
    passed = all(!harmonized_tnk_candidate | stats::complete.cases(.SD)),
    observed = sum(harmonized_tnk_candidate & !stats::complete.cases(.SD)), expected = 0L
  ), by = dataset, .SDcols = paste0("harmonized_score_", program_signatures)],
  all_meta[, .(
    check = "nk_is_distinct_from_t_subtypes",
    passed = all(is.na(harmonized_tcell_subtype) | harmonized_lineage != "NK" | harmonized_tcell_subtype == "NK"),
    observed = sum(harmonized_lineage == "NK" & harmonized_tcell_subtype != "NK", na.rm = TRUE), expected = 0L
  ), by = dataset]
), use.names = TRUE, fill = TRUE)
data.table::fwrite(qc, file.path(run_dir, "qc_tables", "harmonization_validation.tsv"), sep = "\t", na = "NA")
if (any(!qc$passed)) stop("One or more harmonization validation checks failed")

subtype_plot_data <- all_meta[harmonized_lineage == "T_cell", .N, by = .(
  dataset, patient_id_harmonized, harmonized_tcell_subtype
)]
subtype_plot_data[, fraction := N / sum(N), by = .(dataset, patient_id_harmonized)]
p_subtype <- ggplot(subtype_plot_data, aes(patient_id_harmonized, fraction, fill = harmonized_tcell_subtype)) +
  geom_col(width = 0.85) + facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Patient", y = "Fraction within T cells", fill = "Subtype") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(run_dir, "figures", "subtype_composition_by_patient.pdf"), p_subtype, width = 12, height = 5)
ggsave(file.path(run_dir, "figures", "subtype_composition_by_patient.png"), p_subtype, width = 12, height = 5, dpi = 160)

program_plot_data <- data.table::melt(
  patient_summary,
  id.vars = c("dataset", "patient_id"),
  measure.vars = paste0("mean_program_", program_signatures),
  variable.name = "program", value.name = "mean_ucell", variable.factor = FALSE
)
program_plot_data[, program := sub("^mean_program_", "", program)]
p_program <- ggplot(program_plot_data, aes(patient_id, program, fill = mean_ucell)) +
  geom_tile() + facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Patient", y = "Program", fill = "Mean UCell") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(run_dir, "figures", "program_scores_by_patient.pdf"), p_program, width = 12, height = 4.8)
ggsave(file.path(run_dir, "figures", "program_scores_by_patient.png"), p_program, width = 12, height = 4.8, dpi = 160)

gbm_umap <- gbm[harmonized_tnk_candidate %in% TRUE & is.finite(harmonized_tcell_umap_1)]
p_umap <- ggplot(gbm_umap, aes(harmonized_tcell_umap_1, harmonized_tcell_umap_2, color = harmonized_tcell_subtype)) +
  geom_point(size = 0.55, alpha = 0.8) + coord_equal() +
  labs(x = "GBM T/NK UMAP 1", y = "GBM T/NK UMAP 2", color = "Subtype") +
  theme_bw(base_size = 10)
ggsave(file.path(run_dir, "figures", "GBM_tcell_subtype_umap.pdf"), p_umap, width = 7, height = 5.5)
ggsave(file.path(run_dir, "figures", "GBM_tcell_subtype_umap.png"), p_umap, width = 7, height = 5.5, dpi = 180)

writeLines(c(
  paste0("completed_at=", timestamp()),
  paste0("hostname=", Sys.info()[["nodename"]]),
  paste0("validation_checks=", nrow(qc)),
  "validation_status=PASS"
), file.path(run_dir, "provenance", "qc_complete.txt"))
log_message("QC and cross-dataset comparison complete")
