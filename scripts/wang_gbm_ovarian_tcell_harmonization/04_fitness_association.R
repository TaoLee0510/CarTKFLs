#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 04_fitness_association.R CONFIG_YAML RUN_DIR")
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

parse_karyotypes <- function(k, n_chr) {
  pieces <- strsplit(as.character(k), ".", fixed = TRUE)
  valid <- lengths(pieces) == n_chr
  out <- matrix(NA_real_, nrow = length(k), ncol = n_chr)
  if (any(valid)) out[valid, ] <- t(vapply(pieces[valid], as.numeric, numeric(n_chr)))
  colnames(out) <- paste0("chr", seq_len(n_chr))
  out
}

landscape_effects <- function(dataset, row) {
  path <- file.path(
    cfg$fitness$landscape_root, dataset, "ALFAK_fitnessLandscape",
    as.character(row$pm_label), paste0("MINOBS_", row$min_obs),
    as.character(row$PatientID), "landscape.Rds"
  )
  if (!file.exists(path)) {
    return(data.table(dataset = dataset, patient_id = as.character(row$PatientID), status = "landscape_missing", landscape_path = path))
  }
  d <- as.data.table(readRDS(path))
  required <- c("k", "mean", "fq")
  if (length(setdiff(required, colnames(d)))) stop("Landscape lacks required columns: ", path)
  d <- d[fq %in% TRUE & is.finite(mean)]
  minimum_n <- as.integer(cfg$fitness$minimum_frequent_karyotypes)
  if (nrow(d) < minimum_n) {
    return(data.table(dataset = dataset, patient_id = as.character(row$PatientID), status = "too_few_frequent_karyotypes", landscape_path = path))
  }
  copies <- parse_karyotypes(d$k, as.integer(cfg$fitness$chromosomes))
  rows <- lapply(seq_len(ncol(copies)), function(j) {
    x <- copies[, j]
    y <- d$mean
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    n_unique <- length(unique(x))
    if (length(x) < minimum_n || n_unique < as.integer(cfg$fitness$minimum_unique_copy_numbers)) {
      return(data.table(chromosome = j, n_karyotypes = length(x), n_unique_copy_numbers = n_unique, status = "not_estimable"))
    }
    fit <- stats::lm(y ~ x)
    ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
    data.table(
      chromosome = j,
      n_karyotypes = length(x),
      n_unique_copy_numbers = n_unique,
      marginal_slope = unname(stats::coef(fit)[[2L]]),
      standardized_beta = suppressWarnings(stats::cor(x, y, method = "pearson")),
      within_landscape_spearman_rho = unname(ct$estimate),
      within_landscape_spearman_p = ct$p.value,
      status = "estimated"
    )
  })
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, `:=`(
    dataset = dataset,
    patient_id = as.character(row$PatientID),
    pm_label = as.character(row$pm_label),
    pm_value = as.numeric(row$pm_value),
    min_obs = as.integer(row$min_obs),
    selected_fit_correlation = as.numeric(row$correlation),
    selected_fit_r_squared = as.numeric(row$r_squared),
    landscape_path = path
  )]
  out
}

log_message("Deriving chromosome-specific effects from selected patient landscapes")
effect_tables <- list()
selection_audit <- list()
for (dataset in c("GBM", "Ovarian")) {
  best_path <- file.path(cfg$fitness$best_parameter_root, dataset, "best_correlation_results.csv")
  best <- fread(best_path)
  best[, PatientID := as.character(PatientID)]
  if (dataset == "Ovarian") {
    best[, primary_patient_row := !grepl("^EOC153_omentum_", PatientID)]
  } else {
    best[, primary_patient_row := TRUE]
  }
  selection_audit[[dataset]] <- best[, .(dataset = dataset, PatientID, pm_label, pm_value, min_obs, primary_patient_row)]
  primary <- best[primary_patient_row %in% TRUE]
  effect_tables[[dataset]] <- rbindlist(lapply(seq_len(nrow(primary)), function(i) {
    landscape_effects(dataset, primary[i])
  }), use.names = TRUE, fill = TRUE)
}
effects <- rbindlist(effect_tables, use.names = TRUE, fill = TRUE)
fwrite(effects, file.path(run_dir, "association_results", "patient_chromosome_fitness_effects.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(selection_audit), file.path(run_dir, "association_results", "best_parameter_row_audit.tsv"), sep = "\t", na = "NA")

patient <- fread(file.path(run_dir, "sample_summaries", "harmonized_patient_summary.tsv"))
joined <- merge(effects[status == "estimated"], patient, by = c("dataset", "patient_id"), all.x = TRUE)
metric_names <- c("marginal_slope", "standardized_beta", "within_landscape_spearman_rho")
feature_names <- c(
  "total_t_fraction_all_cells",
  grep("^prop_within_t_", colnames(patient), value = TRUE),
  grep("^mean_program_", colnames(patient), value = TRUE)
)
feature_names <- unique(feature_names)

cor_rows <- list()
k <- 0L
for (ds in unique(joined$dataset)) {
  for (chr in sort(unique(joined[dataset == ds]$chromosome))) {
    block <- joined[dataset == ds & chromosome == chr]
    for (metric in metric_names) {
      for (feature in feature_names) {
        keep <- is.finite(block[[metric]]) & is.finite(block[[feature]])
        k <- k + 1L
        if (sum(keep) < 5L || length(unique(block[[feature]][keep])) < 2L) {
          cor_rows[[k]] <- data.table(dataset = ds, chromosome = chr, chromosome_effect_metric = metric, feature = feature, n_patients = sum(keep), spearman_rho = NA_real_, p_value = NA_real_)
        } else {
          ct <- suppressWarnings(stats::cor.test(block[[metric]][keep], block[[feature]][keep], method = "spearman", exact = FALSE))
          cor_rows[[k]] <- data.table(dataset = ds, chromosome = chr, chromosome_effect_metric = metric, feature = feature, n_patients = sum(keep), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
        }
      }
    }
  }
}
associations <- rbindlist(cor_rows)
associations[, fdr_bh := stats::p.adjust(p_value, method = "BH"), by = .(dataset, chromosome_effect_metric)]
fwrite(associations, file.path(run_dir, "association_results", "chromosome_feature_spearman.tsv"), sep = "\t", na = "NA")
primary <- associations[feature == "total_t_fraction_all_cells"]
fwrite(primary, file.path(run_dir, "association_results", "primary_total_t_associations.tsv"), sep = "\t", na = "NA")

adjust_rows <- list()
k <- 0L
composition_state_features <- setdiff(feature_names, "total_t_fraction_all_cells")
for (ds in unique(joined$dataset)) {
  for (chr in sort(unique(joined[dataset == ds]$chromosome))) {
    block <- joined[dataset == ds & chromosome == chr]
    for (metric in metric_names) {
      for (feature in composition_state_features) {
        keep <- stats::complete.cases(block[, c(metric, "total_t_fraction_all_cells", feature), with = FALSE])
        if (sum(keep) < 7L || length(unique(block[[feature]][keep])) < 2L) next
        model_data <- data.frame(
          y = as.numeric(scale(block[[metric]][keep])),
          total_t = as.numeric(scale(block$total_t_fraction_all_cells[keep])),
          added_feature = as.numeric(scale(block[[feature]][keep]))
        )
        if (!all(stats::complete.cases(model_data))) next
        base_fit <- stats::lm(y ~ total_t, data = model_data)
        adjusted_fit <- stats::lm(y ~ total_t + added_feature, data = model_data)
        base_coef <- summary(base_fit)$coefficients
        adjusted_coef <- summary(adjusted_fit)$coefficients
        k <- k + 1L
        adjust_rows[[k]] <- data.table(
          dataset = ds, chromosome = chr, chromosome_effect_metric = metric,
          added_feature = feature, n_patients = nrow(model_data),
          base_total_t_beta = base_coef["total_t", "Estimate"],
          base_total_t_p = base_coef["total_t", "Pr(>|t|)"],
          adjusted_total_t_beta = adjusted_coef["total_t", "Estimate"],
          adjusted_total_t_p = adjusted_coef["total_t", "Pr(>|t|)"],
          added_feature_beta = adjusted_coef["added_feature", "Estimate"],
          added_feature_p = adjusted_coef["added_feature", "Pr(>|t|)"],
          base_r_squared = summary(base_fit)$r.squared,
          adjusted_r_squared = summary(adjusted_fit)$r.squared,
          delta_r_squared = summary(adjusted_fit)$r.squared - summary(base_fit)$r.squared,
          total_t_beta_attenuation_fraction = 1 - abs(adjusted_coef["total_t", "Estimate"]) / max(abs(base_coef["total_t", "Estimate"]), .Machine$double.eps)
        )
      }
    }
  }
}
adjustments <- if (length(adjust_rows)) rbindlist(adjust_rows) else data.table()
if (nrow(adjustments)) {
  adjustments[, added_feature_fdr_bh := stats::p.adjust(added_feature_p, method = "BH"),
              by = .(dataset, chromosome_effect_metric)]
  adjustments[, adjusted_total_t_fdr_bh := stats::p.adjust(adjusted_total_t_p, method = "BH"),
              by = .(dataset, chromosome_effect_metric)]
}
fwrite(adjustments, file.path(run_dir, "association_results", "one_feature_adjustment_models.tsv"), sep = "\t", na = "NA")

common <- dcast(
  associations[is.finite(spearman_rho)],
  chromosome + chromosome_effect_metric + feature ~ dataset,
  value.var = "spearman_rho"
)
if (all(c("GBM", "Ovarian") %in% colnames(common))) {
  common[, `:=`(
    opposite_direction = is.finite(GBM) & is.finite(Ovarian) & sign(GBM) != sign(Ovarian),
    absolute_direction_difference = abs(GBM - Ovarian)
  )]
  setorder(common, -opposite_direction, -absolute_direction_difference)
}
fwrite(common, file.path(run_dir, "association_results", "cross_cancer_direction_comparison.tsv"), sep = "\t", na = "NA")

plot_data <- associations[chromosome_effect_metric == "marginal_slope" & is.finite(spearman_rho)]
plot_data[, feature_label := sub("^mean_program_", "", sub("^prop_within_t_", "prop_", feature))]
p <- ggplot(plot_data, aes(factor(chromosome), feature_label, fill = spearman_rho)) +
  geom_tile() + facet_grid(. ~ dataset) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  labs(x = "Chromosome", y = "T-cell composition/state feature", fill = "Spearman rho") +
  theme_bw(base_size = 9) + theme(axis.text.y = element_text(size = 7))
ggsave(file.path(run_dir, "figures", "chromosome_fitness_feature_associations.pdf"), p, width = 13, height = 7)
ggsave(file.path(run_dir, "figures", "chromosome_fitness_feature_associations.png"), p, width = 13, height = 7, dpi = 170)

writeLines(c(
  paste0("completed_at=", timestamp()),
  paste0("hostname=", Sys.info()[["nodename"]]),
  paste0("patient_chromosome_rows=", nrow(effects)),
  paste0("feature_association_rows=", nrow(associations)),
  "interpretation=cross_patient_exploratory_association_not_causal"
), file.path(run_dir, "provenance", "fitness_association_complete.txt"))
log_message("Fitness association analysis complete")
