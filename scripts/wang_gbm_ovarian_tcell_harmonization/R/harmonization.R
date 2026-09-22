split_genes <- function(x) {
  y <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  unique(toupper(y[nzchar(y)]))
}

read_marker_programs <- function(path) {
  x <- read_tsv(path)
  required <- c("signature", "level", "label", "genes", "citation", "pmid", "doi", "url")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) stop("Marker table lacks: ", paste(missing, collapse = ", "))
  if (anyDuplicated(x$signature)) stop("Marker signatures must be unique")
  x
}

select_expression_assay <- function(obj, preferences) {
  available <- Seurat::Assays(obj)
  for (assay in preferences) {
    if (!assay %in% available) next
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, slot = "data"),
      error = function(e) NULL
    )
    if (!is.null(mat) && nrow(mat) > 0L && ncol(mat) > 0L && Matrix::nnzero(mat) > 0L) {
      return(assay)
    }
  }
  stop("No non-empty normalized expression assay found among: ", paste(preferences, collapse = ", "))
}

matched_signature_sets <- function(marker_table, object_genes) {
  object_upper <- toupper(object_genes)
  keep_first <- !duplicated(object_upper)
  gene_lookup <- setNames(object_genes[keep_first], object_upper[keep_first])
  out <- lapply(marker_table$genes, function(g) {
    requested <- split_genes(g)
    unname(gene_lookup[intersect(requested, names(gene_lookup))])
  })
  names(out) <- marker_table$signature
  out
}

signature_coverage_table <- function(marker_table, matched_sets, dataset) {
  requested <- lapply(marker_table$genes, split_genes)
  data.frame(
    dataset = dataset,
    signature = marker_table$signature,
    level = marker_table$level,
    label = marker_table$label,
    requested_genes = vapply(requested, paste, collapse = ";", character(1)),
    observed_genes = vapply(matched_sets, function(x) paste(toupper(x), collapse = ";"), character(1)),
    n_requested = lengths(requested),
    n_observed = lengths(matched_sets),
    coverage_fraction = lengths(matched_sets) / pmax(1L, lengths(requested)),
    citation = marker_table$citation,
    pmid = marker_table$pmid,
    doi = marker_table$doi,
    url = marker_table$url,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

select_tnk_cells <- function(obj, dataset) {
  meta <- obj@meta.data
  if (identical(dataset, "GBM")) {
    if (!"Normal_celltype" %in% colnames(meta)) stop("GBM lacks Normal_celltype")
    keep <- safe_character(meta$Normal_celltype) == "T cell"
  } else if (identical(dataset, "Ovarian")) {
    required <- c("cell_type_major_final", "cell_type_subtype_final")
    if (length(setdiff(required, colnames(meta)))) stop("Ovarian lacks final hierarchical annotation")
    major <- safe_character(meta$cell_type_major_final)
    broad <- if ("cell_type_broad" %in% colnames(meta)) safe_character(meta$cell_type_broad) else ""
    keep <- major %in% c("T_cell", "NK_cell") |
      grepl("(^|[:| ])T_cell($|[| ])|(^|[:| ])NK_cell($|[| ])", broad)
    if ("annotation_analysis_disposition_final" %in% colnames(meta)) {
      keep <- keep & safe_character(meta$annotation_analysis_disposition_final) != "exclude_from_biological_analysis"
    }
  } else {
    stop("Unsupported dataset: ", dataset)
  }
  rownames(meta)[which(keep)]
}

score_tnk_cells <- function(obj, cells, marker_table, cfg, dataset) {
  assay <- select_expression_assay(obj, unlist(cfg$annotation$assay_preference))
  expr <- Seurat::GetAssayData(obj, assay = assay, slot = "data")[, cells, drop = FALSE]
  matched <- matched_signature_sets(marker_table, rownames(expr))
  coverage <- signature_coverage_table(marker_table, matched, dataset)
  min_fraction <- as.numeric(cfg$annotation$minimum_signature_coverage_fraction)
  min_genes <- as.integer(cfg$annotation$minimum_signature_genes)
  failed <- coverage$coverage_fraction < min_fraction | coverage$n_observed < min_genes
  if (any(failed)) {
    stop(
      dataset, " has insufficient gene coverage for signatures: ",
      paste(coverage$signature[failed], collapse = ", ")
    )
  }
  log_message("Scoring ", length(cells), " ", dataset, " T/NK candidate cells with UCell")
  scores <- UCell::ScoreSignatures_UCell(
    matrix = expr,
    features = matched,
    maxRank = as.integer(cfg$runtime$ucell_max_rank),
    ncores = as.integer(cfg$runtime$workers),
    force.gc = TRUE
  )
  scores <- as.data.frame(scores, check.names = FALSE)
  colnames(scores) <- sub("_UCell$", "", colnames(scores))
  if (!identical(rownames(scores), cells)) scores <- scores[cells, , drop = FALSE]
  list(scores = scores, coverage = coverage, assay = assay)
}

ovarian_harmonized_labels <- function(meta, cells) {
  subtype <- safe_character(meta[cells, "cell_type_subtype_final"])
  map <- c(
    CD8_T_cell = "CD8",
    CD4_T_cell = "Conventional_CD4",
    Regulatory_T_cell = "Treg",
    T_cell_unresolved = "Unresolved_T",
    NK_cell = "NK",
    Immune_unresolved = "Unresolved_T",
    Unresolved = "Unresolved_T"
  )
  label <- unname(map[subtype])
  label[is.na(label)] <- "Unresolved_T"
  confidence <- if ("annotation_confidence_final" %in% colnames(meta)) {
    safe_character(meta[cells, "annotation_confidence_final"])
  } else {
    rep("not_recorded", length(cells))
  }
  reason <- if ("annotation_qc_flag_final" %in% colnames(meta)) {
    safe_character(meta[cells, "annotation_qc_flag_final"])
  } else {
    rep("approved_ovarian_mapping", length(cells))
  }
  data.frame(
    cell_id = cells,
    harmonized_tcell_subtype = label,
    harmonized_subtype_confidence = confidence,
    harmonized_uncertainty_reason = reason,
    harmonized_assignment_source = "approved_ovarian_hierarchical_v4_mapping",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

fit_ovarian_reference <- function(scores, labels) {
  feature_names <- c("T_lineage", "NK_lineage", "CD8_subtype", "CD4_conventional_subtype", "Treg_subtype")
  classes <- c("CD8", "Conventional_CD4", "Treg", "NK")
  idx <- labels %in% classes & stats::complete.cases(scores[, feature_names, drop = FALSE])
  if (sum(idx) < 100L) stop("Too few ovarian reference cells")
  x <- as.matrix(scores[idx, feature_names, drop = FALSE])
  y <- labels[idx]
  center <- apply(x, 2L, stats::median, na.rm = TRUE)
  scale <- apply(x, 2L, stats::mad, na.rm = TRUE)
  fallback <- apply(x, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- fallback[!is.finite(scale) | scale <= 0]
  scale[!is.finite(scale) | scale <= 0] <- 1
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  centroids <- do.call(rbind, lapply(classes, function(cl) {
    colMeans(z[y == cl, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(centroids) <- classes
  own_dist <- numeric(nrow(z))
  for (cl in classes) {
    ii <- which(y == cl)
    own_dist[ii] <- sqrt(rowSums((z[ii, , drop = FALSE] - matrix(
      centroids[cl, ], nrow = length(ii), ncol = ncol(z), byrow = TRUE
    ))^2))
  }
  radii <- vapply(classes, function(cl) {
    as.numeric(stats::quantile(own_dist[y == cl], 0.95, na.rm = TRUE, names = FALSE))
  }, numeric(1))
  radii[!is.finite(radii) | radii <= 0] <- 0.5
  list(
    feature_names = feature_names,
    center = center,
    scale = scale,
    centroids = centroids,
    radii = radii,
    n_reference = table(factor(y, levels = classes))
  )
}

reference_to_table <- function(reference) {
  out <- do.call(rbind, lapply(rownames(reference$centroids), function(cl) {
    base <- data.frame(
      class = cl,
      n_reference = as.integer(reference$n_reference[[cl]]),
      radius_95pct = unname(reference$radii[[cl]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    cbind(base, as.data.frame(as.list(reference$centroids[cl, ]), check.names = FALSE))
  }))
  rownames(out) <- NULL
  out
}

recluster_gbm_tcells <- function(obj, cells, cfg) {
  t_obj <- obj[, cells]
  Seurat::DefaultAssay(t_obj) <- "RNA"
  set.seed(as.integer(cfg$seed))
  t_obj <- Seurat::NormalizeData(t_obj, verbose = FALSE)
  t_obj <- Seurat::FindVariableFeatures(t_obj, nfeatures = min(2000L, nrow(t_obj)), verbose = FALSE)
  t_obj <- Seurat::ScaleData(t_obj, features = Seurat::VariableFeatures(t_obj), verbose = FALSE)
  npcs <- min(30L, length(Seurat::VariableFeatures(t_obj)) - 1L, ncol(t_obj) - 1L)
  t_obj <- Seurat::RunPCA(t_obj, npcs = npcs, verbose = FALSE)
  dims_use <- seq_len(min(as.integer(cfg$annotation$gbm_recluster_dims), npcs))
  reduction <- "pca"
  if ("PatientID" %in% colnames(t_obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
    t_obj <- harmony::RunHarmony(
      object = t_obj, group.by.vars = "PatientID", reduction.use = "pca",
      dims.use = dims_use, plot_convergence = FALSE, verbose = FALSE
    )
    reduction <- "harmony"
  }
  t_obj <- Seurat::FindNeighbors(t_obj, reduction = reduction, dims = dims_use, verbose = FALSE)
  t_obj <- Seurat::FindClusters(
    t_obj, resolution = as.numeric(cfg$annotation$gbm_recluster_resolution),
    random.seed = as.integer(cfg$seed), verbose = FALSE
  )
  t_obj <- Seurat::RunUMAP(
    t_obj, reduction = reduction, dims = dims_use, seed.use = as.integer(cfg$seed), verbose = FALSE
  )
  t_obj
}

classify_gbm_clusters <- function(scores, clusters, reference, cfg) {
  feature_names <- reference$feature_names
  cluster_ids <- sort(unique(as.character(clusters)))
  rows <- lapply(cluster_ids, function(cl) {
    cells <- names(clusters)[as.character(clusters) == cl]
    med <- vapply(feature_names, function(f) stats::median(scores[cells, f], na.rm = TRUE), numeric(1))
    z <- (med - reference$center) / reference$scale
    distances <- sqrt(rowSums((reference$centroids - matrix(
      z, nrow = nrow(reference$centroids), ncol = length(z), byrow = TRUE
    ))^2))
    ord <- order(distances)
    best <- names(distances)[ord[[1L]]]
    d1 <- unname(distances[ord[[1L]]])
    d2 <- unname(distances[ord[[2L]]])
    rel_margin <- (d2 - d1) / max(d2, .Machine$double.eps)
    within_radius <- d1 <= reference$radii[[best]] * as.numeric(cfg$annotation$centroid_radius_multiplier)
    lineage_ok <- if (identical(best, "NK")) {
      med[["NK_lineage"]] > med[["T_lineage"]]
    } else {
      med[["T_lineage"]] >= med[["NK_lineage"]] - as.numeric(cfg$annotation$score_margin_moderate)
    }
    enough_cells <- length(cells) >= as.integer(cfg$annotation$minimum_cluster_cells)
    accepted <- within_radius && lineage_ok && enough_cells &&
      rel_margin >= as.numeric(cfg$annotation$centroid_relative_margin_moderate)
    label <- if (accepted) best else "Unresolved_T"
    confidence <- if (!accepted) {
      "low"
    } else if (rel_margin >= as.numeric(cfg$annotation$centroid_relative_margin_high) &&
               d1 <= reference$radii[[best]]) {
      "high"
    } else {
      "moderate"
    }
    reason <- paste(
      c(
        if (!within_radius) "outside_ovarian_reference_radius" else NULL,
        if (!lineage_ok) "lineage_score_conflict" else NULL,
        if (!enough_cells) "cluster_below_minimum_cells" else NULL,
        if (rel_margin < as.numeric(cfg$annotation$centroid_relative_margin_moderate)) "low_centroid_margin" else NULL
      ),
      collapse = ";"
    )
    if (!nzchar(reason)) reason <- "passed_reference_centroid_rules"
    base <- data.frame(
      harmonized_cluster = cl,
      n_cells = length(cells),
      nearest_reference_class = best,
      nearest_distance = d1,
      second_distance = d2,
      relative_distance_margin = rel_margin,
      reference_radius_95pct = unname(reference$radii[[best]]),
      within_reference_radius = within_radius,
      lineage_check_pass = lineage_ok,
      harmonized_tcell_subtype = label,
      harmonized_subtype_confidence = confidence,
      harmonized_uncertainty_reason = reason,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    cbind(base, as.data.frame(as.list(med), check.names = FALSE))
  })
  do.call(rbind, rows)
}

derive_program_states <- function(scores, marker_table, cfg) {
  program_names <- marker_table$signature[marker_table$level == "program"]
  labels <- setNames(marker_table$label, marker_table$signature)
  raw <- as.matrix(scores[, program_names, drop = FALSE])
  z <- raw
  for (j in seq_len(ncol(raw))) {
    center <- stats::median(raw[, j], na.rm = TRUE)
    spread <- stats::mad(raw[, j], na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) spread <- stats::sd(raw[, j], na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) spread <- 1
    z[, j] <- (raw[, j] - center) / spread
  }
  colnames(z) <- paste0(program_names, "_robust_z")
  threshold <- as.numeric(cfg$annotation$state_z_support)
  margin_threshold <- as.numeric(cfg$annotation$state_z_margin)
  supported <- apply(z, 1L, function(v) {
    hit <- which(is.finite(v) & v >= threshold)
    if (!length(hit)) return("None_supported")
    paste(unname(labels[sub("_robust_z$", "", names(v)[hit])]), collapse = ";")
  })
  primary <- apply(z, 1L, function(v) {
    ord <- order(v, decreasing = TRUE, na.last = NA)
    if (!length(ord) || v[ord[[1L]]] < threshold) return("None_supported")
    if (length(ord) > 1L && (v[ord[[1L]]] - v[ord[[2L]]]) < margin_threshold) {
      return("Mixed_supported_programs")
    }
    unname(labels[sub("_robust_z$", "", names(v)[ord[[1L]]])])
  })
  list(z = as.data.frame(z, check.names = FALSE), supported = supported, primary = primary)
}

add_harmonized_metadata <- function(obj, cells, scores, assignments, states, cluster_ids, dataset) {
  meta <- obj@meta.data
  all_cells <- rownames(meta)
  add_character <- function(values, default = NA_character_) {
    out <- rep(default, length(all_cells)); names(out) <- all_cells
    out[cells] <- as.character(values)
    out
  }
  add_numeric <- function(values) {
    out <- rep(NA_real_, length(all_cells)); names(out) <- all_cells
    out[cells] <- as.numeric(values)
    out
  }
  meta$harmonized_tnk_candidate <- all_cells %in% cells
  meta$harmonized_lineage <- add_character(ifelse(assignments$harmonized_tcell_subtype == "NK", "NK", "T_cell"))
  meta$harmonized_tcell_subtype <- add_character(assignments$harmonized_tcell_subtype)
  meta$harmonized_subtype_confidence <- add_character(assignments$harmonized_subtype_confidence)
  meta$harmonized_uncertainty_reason <- add_character(assignments$harmonized_uncertainty_reason)
  meta$harmonized_assignment_source <- add_character(assignments$harmonized_assignment_source)
  meta$harmonized_cluster <- add_character(cluster_ids)
  meta$harmonized_state_primary <- add_character(states$primary)
  meta$harmonized_state_supported <- add_character(states$supported)
  for (sig in colnames(scores)) {
    meta[[paste0("harmonized_score_", sig)]] <- add_numeric(scores[[sig]])
  }
  for (sig in colnames(states$z)) {
    meta[[paste0("harmonized_score_", sig)]] <- add_numeric(states$z[[sig]])
  }
  meta$harmonized_annotation_dataset <- dataset
  obj@meta.data <- meta
  obj
}
