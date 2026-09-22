#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: test_helpers.R SCRIPT_DIR")
script_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source(file.path(script_dir, "R", "common.R"))
source(file.path(script_dir, "R", "harmonization.R"))

set.seed(1)
score_names <- c(
  "T_lineage", "NK_lineage", "CD8_subtype", "CD4_conventional_subtype", "Treg_subtype",
  "cytotoxicity", "activation", "naive_memory", "exhaustion_dysfunction", "interferon_response", "stress"
)
scores <- as.data.frame(matrix(runif(400L * length(score_names)), nrow = 400L))
colnames(scores) <- score_names
rownames(scores) <- paste0("cell", seq_len(nrow(scores)))
labels <- rep(c("CD8", "Conventional_CD4", "Treg", "NK"), each = 100L)
scores[labels == "CD8", "CD8_subtype"] <- scores[labels == "CD8", "CD8_subtype"] + 2
scores[labels == "Conventional_CD4", "CD4_conventional_subtype"] <- scores[labels == "Conventional_CD4", "CD4_conventional_subtype"] + 2
scores[labels == "Treg", "Treg_subtype"] <- scores[labels == "Treg", "Treg_subtype"] + 2
scores[labels == "NK", "NK_lineage"] <- scores[labels == "NK", "NK_lineage"] + 2
ref <- fit_ovarian_reference(scores, labels)
stopifnot(nrow(reference_to_table(ref)) == 4L, all(is.finite(ref$radii)), all(ref$radii > 0))

markers <- read_marker_programs(file.path(script_dir, "marker_programs.tsv"))
states <- derive_program_states(scores, markers, list(annotation = list(state_z_support = 1, state_z_margin = 0.5)))
stopifnot(nrow(states$z) == nrow(scores), length(states$primary) == nrow(scores))

suppressPackageStartupMessages(library(Seurat))
toy_counts <- matrix(rpois(60, lambda = 2), nrow = 6, dimnames = list(
  paste0("gene", seq_len(6)), paste0("toy", seq_len(10))
))
toy <- CreateSeuratObject(toy_counts)
toy$Normal_celltype <- c(rep("T cell", 3), rep("Other", 7))
selected <- select_tnk_cells(toy, "GBM")
stopifnot(identical(selected, colnames(toy)[seq_len(3)]))
toy_subset <- toy[, selected]
stopifnot(ncol(toy_subset) == 3L)
cat("helper_tests=PASS\n")
