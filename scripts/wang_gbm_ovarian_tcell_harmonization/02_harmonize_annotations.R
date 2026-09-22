#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 02_harmonize_annotations.R CONFIG_YAML RUN_DIR")
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
run_dir <- normalizePath(args[[2L]], mustWork = TRUE)
script_dir <- dirname(config_path)
source(file.path(script_dir, "R", "common.R"))
source(file.path(script_dir, "R", "harmonization.R"))

cfg <- read_analysis_config(config_path)
ensure_run_directories(run_dir)
activate_task_library(run_dir, cfg)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(data.table)
})

if (!identical(Sys.info()[["nodename"]], cfg$required_hostname)) {
  stop("Required host is ", cfg$required_hostname, "; observed ", Sys.info()[["nodename"]])
}
set.seed(as.integer(cfg$seed))
markers <- read_marker_programs(file.path(script_dir, "marker_programs.tsv"))

harmonization_checkpoint <- c(
  sentinel = file.path(run_dir, "provenance", "harmonization_complete.txt"),
  ovarian_object = file.path(run_dir, "annotated_objects", "Ovarian_harmonized_tcell_annotations.Rds"),
  gbm_object = file.path(run_dir, "annotated_objects", "GBM_harmonized_tcell_annotations.Rds"),
  ovarian_metadata = file.path(run_dir, "cell_metadata", "Ovarian_harmonized_cell_metadata.tsv.gz"),
  gbm_metadata = file.path(run_dir, "cell_metadata", "GBM_harmonized_cell_metadata.tsv.gz"),
  reference = file.path(run_dir, "provenance", "ovarian_subtype_reference.rds")
)
if (all(file.exists(harmonization_checkpoint)) &&
    all(file.info(harmonization_checkpoint[c("ovarian_object", "gbm_object")])$size > 1e6)) {
  log_message("Reusing completed GBM and ovarian harmonization checkpoint")
  quit(save = "no", status = 0L)
}

ovarian_checkpoint <- c(
  object = file.path(run_dir, "annotated_objects", "Ovarian_harmonized_tcell_annotations.Rds"),
  metadata = file.path(run_dir, "cell_metadata", "Ovarian_harmonized_cell_metadata.tsv.gz"),
  reference = file.path(run_dir, "provenance", "ovarian_subtype_reference.rds"),
  coverage = file.path(run_dir, "qc_tables", "Ovarian_signature_gene_coverage.tsv")
)
if (all(file.exists(ovarian_checkpoint)) && file.info(ovarian_checkpoint[["object"]])$size > 1e6) {
  log_message("Reusing completed ovarian harmonization checkpoint")
  reference <- readRDS(ovarian_checkpoint[["reference"]])
  ovarian_candidate_flag <- data.table::fread(
    ovarian_checkpoint[["metadata"]], select = "harmonized_tnk_candidate"
  )[["harmonized_tnk_candidate"]]
  ovarian_cells_n <- sum(ovarian_candidate_flag %in% TRUE, na.rm = TRUE)
  rm(ovarian_candidate_flag)
} else {
  log_message("Reading ovarian reference object")
  ovarian <- readRDS(cfg$inputs$Ovarian)
  ovarian_cells <- select_tnk_cells(ovarian, "Ovarian")
  ovarian_scored <- score_tnk_cells(ovarian, ovarian_cells, markers, cfg, "Ovarian")
  ovarian_assign <- ovarian_harmonized_labels(ovarian@meta.data, ovarian_cells)
  rownames(ovarian_assign) <- ovarian_assign$cell_id
  ovarian_assign <- ovarian_assign[ovarian_cells, , drop = FALSE]
  ovarian_states <- derive_program_states(ovarian_scored$scores, markers, cfg)
  ovarian_clusters <- safe_character(ovarian@meta.data[ovarian_cells, "annotation_cluster_final"])
  names(ovarian_clusters) <- ovarian_cells

  reference <- fit_ovarian_reference(
    ovarian_scored$scores,
    ovarian_assign$harmonized_tcell_subtype
  )
  write_tsv(reference_to_table(reference), file.path(run_dir, "qc_tables", "ovarian_reference_centroids.tsv"))
  saveRDS(reference, ovarian_checkpoint[["reference"]])
  write_tsv(ovarian_scored$coverage, ovarian_checkpoint[["coverage"]])

  ovarian <- add_harmonized_metadata(
    ovarian, ovarian_cells, ovarian_scored$scores, ovarian_assign,
    ovarian_states, ovarian_clusters, "Ovarian"
  )
  ovarian_meta <- data.table::as.data.table(ovarian@meta.data, keep.rownames = "cell_id")
  ovarian_meta[, dataset := "Ovarian"]
  data.table::fwrite(
    ovarian_meta, ovarian_checkpoint[["metadata"]],
    sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
  )
  log_message("Saving harmonized ovarian object")
  saveRDS(ovarian, ovarian_checkpoint[["object"]], compress = TRUE)
  ovarian_cells_n <- length(ovarian_cells)
  rm(ovarian, ovarian_meta, ovarian_states)
  invisible(gc(verbose = FALSE))
}

log_message("Reading GBM object")
gbm <- readRDS(cfg$inputs$GBM)
gbm_cells <- select_tnk_cells(gbm, "GBM")
gbm_scored <- score_tnk_cells(gbm, gbm_cells, markers, cfg, "GBM")
write_tsv(gbm_scored$coverage, file.path(run_dir, "qc_tables", "GBM_signature_gene_coverage.tsv"))

log_message("Running GBM T-cell-specific reclustering")
gbm_t <- recluster_gbm_tcells(gbm, gbm_cells, cfg)
gbm_clusters <- safe_character(gbm_t$seurat_clusters)
names(gbm_clusters) <- colnames(gbm_t)
gbm_clusters <- gbm_clusters[gbm_cells]
gbm_scored$scores <- gbm_scored$scores[gbm_cells, , drop = FALSE]
cluster_calls <- classify_gbm_clusters(gbm_scored$scores, gbm_clusters, reference, cfg)
write_tsv(cluster_calls, file.path(run_dir, "qc_tables", "GBM_cluster_subtype_assignments.tsv"))
call_index <- match(gbm_clusters, cluster_calls$harmonized_cluster)
gbm_assign <- data.frame(
  cell_id = gbm_cells,
  harmonized_tcell_subtype = cluster_calls$harmonized_tcell_subtype[call_index],
  harmonized_subtype_confidence = cluster_calls$harmonized_subtype_confidence[call_index],
  harmonized_uncertainty_reason = cluster_calls$harmonized_uncertainty_reason[call_index],
  harmonized_assignment_source = "gbm_tcell_harmony_cluster_to_ovarian_reference_centroid",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gbm_states <- derive_program_states(gbm_scored$scores, markers, cfg)
gbm <- add_harmonized_metadata(
  gbm, gbm_cells, gbm_scored$scores, gbm_assign,
  gbm_states, gbm_clusters, "GBM"
)
umap <- Seurat::Embeddings(gbm_t, "umap")
gbm@meta.data$harmonized_tcell_umap_1 <- NA_real_
gbm@meta.data$harmonized_tcell_umap_2 <- NA_real_
gbm@meta.data[rownames(umap), "harmonized_tcell_umap_1"] <- umap[, 1L]
gbm@meta.data[rownames(umap), "harmonized_tcell_umap_2"] <- umap[, 2L]
saveRDS(gbm_t, file.path(run_dir, "audit", "GBM_tcell_reclustered_subset.Rds"), compress = TRUE)

gbm_meta <- data.table::as.data.table(gbm@meta.data, keep.rownames = "cell_id")
gbm_meta[, dataset := "GBM"]
data.table::fwrite(
  gbm_meta,
  file.path(run_dir, "cell_metadata", "GBM_harmonized_cell_metadata.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)
log_message("Saving harmonized GBM object")
saveRDS(
  gbm,
  file.path(run_dir, "annotated_objects", "GBM_harmonized_tcell_annotations.Rds"),
  compress = TRUE
)

writeLines(
  c(
    paste0("completed_at=", timestamp()),
    paste0("hostname=", Sys.info()[["nodename"]]),
    paste0("R_version=", R.version.string),
    paste0("Seurat_version=", as.character(utils::packageVersion("Seurat"))),
    paste0("UCell_version=", as.character(utils::packageVersion("UCell"))),
    paste0("Ovarian_TNK_cells=", ovarian_cells_n),
    paste0("GBM_TNK_cells=", length(gbm_cells))
  ),
  file.path(run_dir, "provenance", "harmonization_complete.txt")
)
log_message("Harmonization complete")
