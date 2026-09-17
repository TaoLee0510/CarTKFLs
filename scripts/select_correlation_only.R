args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: select_correlation_only.R CONFIG SOURCE_RUN_DIR VARIANT_RUN_DIR")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))

cfg <- read_config(args[[1]])
source_dir <- normalizePath(args[[2]], mustWork = TRUE)
run_dir <- normalizePath(args[[3]], mustWork = TRUE)
if (identical(source_dir, run_dir) ||
    !identical(dirname(run_dir), source_dir) ||
    !identical(basename(run_dir), "correlation_only")) {
  stop("The correlation-only run must be SOURCE_RUN_DIR/correlation_only")
}
if (file.exists(file.path(run_dir, "manifests", "selected_samples.tsv"))) {
  stop("Correlation-only selection already exists; refusing to overwrite it")
}
if (!file.exists(file.path(source_dir, "fit_audit.tsv"))) {
  stop("Source fit audit is missing")
}

link_source <- function(source, destination) {
  if (!file.exists(source)) stop("Missing source: ", source)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination) || nzchar(Sys.readlink(destination))) {
    if (!identical(normalizePath(destination), normalizePath(source))) {
      stop("Existing link points elsewhere: ", destination)
    }
  } else if (!file.symlink(source, destination)) {
    stop("Cannot link ", source, " to ", destination)
  }
}
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(value, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Cannot install ", path)
}
atomic_copy <- function(source, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  if (!file.copy(source, tmp, overwrite = TRUE) || !file.rename(tmp, path)) {
    stop("Cannot copy ", source, " to ", path)
  }
}

source_cfg <- yaml::read_yaml(file.path(source_dir, "downstream.yaml"))
analysis_id <- as.character(source_cfg$workflow$analysis_id)
if (length(analysis_id) != 1L || !nzchar(analysis_id)) {
  stop("Source downstream analysis_id is missing")
}
source_results <- file.path(source_dir, "analysis", "results", "OptimizedParameters")
original_path <- file.path(source_dir, "manifests", "selected_samples.tsv")
original <- utils::read.delim(original_path, stringsAsFactors = FALSE, check.names = FALSE)
required_manifest <- c("cancer_type", "PatientID", "pm_label", "min_obs")
if (!all(required_manifest %in% names(original))) {
  stop("Original selection manifest is malformed")
}
expected <- expand.grid(
  PatientID = names(cfg$patients),
  high_cn = as.integer(unlist(cfg$high_cn)),
  stringsAsFactors = FALSE
)
expected$cancer_type <- vapply(expected$high_cn, cancer_type, character(1))
expected_key <- paste(expected$cancer_type, expected$PatientID, sep = "\r")
original_key <- paste(original$cancer_type, original$PatientID, sep = "\r")
if (anyDuplicated(original_key) || !all(original_key %in% expected_key)) {
  stop("Original selection has duplicated or unexpected patient/mapping pairs")
}
if (nrow(original) != 5L) {
  stop("Original selection count changed; review SELECTION_METHOD.md")
}

link_source(file.path(source_dir, "analysis", "data"),
            file.path(run_dir, "analysis", "data"))
link_source(file.path(source_dir, "sample_results", "kfl"),
            file.path(run_dir, "sample_results", "kfl"))
link_source(file.path(source_dir, "fit_audit.tsv"),
            file.path(run_dir, "fit_audit.tsv"))

variant_cfg <- source_cfg
variant_cfg$Data_path <- file.path(run_dir, "analysis")
variant_cfg$workflow$analysis_id <- paste0(analysis_id, "_PearsonOnly")
variant_cfg$workflow$cart_selection_method <- "max_finite_pearson_correlation"
variant_cfg$workflow$cart_selection_source_run <- source_dir
config_path <- file.path(run_dir, "downstream.yaml")
tmp_cfg <- tempfile(pattern = ".downstream_", tmpdir = run_dir)
yaml::write_yaml(variant_cfg, tmp_cfg)
if (!file.rename(tmp_cfg, config_path)) {
  unlink(tmp_cfg)
  stop("Cannot install ", config_path)
}

selected_by_ct <- list()
selected_rows <- list()
comparison_rows <- list()
total_evaluated <- 0L
for (ct in unique(expected$cancer_type)) {
  source_metrics <- file.path(source_results, ct, "KFL_metrics.Rds")
  metrics <- readRDS(source_metrics)
  all_results <- metrics$all_results
  required_metrics <- c(
    "PatientID", "pm_label", "pm_value", "min_obs", "correlation",
    "sig_flag", "total_points", "p_value", "huber_perm_p",
    "huber_slope", "signed_r2_Huber_weighted"
  )
  if (!is.data.frame(all_results) ||
      !all(required_metrics %in% names(all_results))) {
    stop("Incomplete KFL metrics: ", source_metrics)
  }
  ids <- expected$PatientID[expected$cancer_type == ct]
  if (!setequal(unique(as.character(all_results$PatientID)), ids)) {
    stop("KFL metrics patient coverage mismatch: ", source_metrics)
  }
  if (anyDuplicated(paste(all_results$PatientID, all_results$pm_label,
                          all_results$min_obs, sep = "\r"))) {
    stop("Duplicated PM/MINOBS rows: ", source_metrics)
  }
  total_evaluated <- total_evaluated + nrow(all_results)

  picks <- lapply(ids, function(pid) {
    group <- all_results[as.character(all_results$PatientID) == pid, , drop = FALSE]
    correlation <- suppressWarnings(as.numeric(group$correlation))
    pm <- suppressWarnings(as.numeric(group$pm_value))
    min_obs <- suppressWarnings(as.integer(group$min_obs))
    eligible <- which(is.finite(correlation) & is.finite(pm) & !is.na(min_obs))
    if (!length(eligible)) stop("No finite Pearson correlation for ", ct, "/", pid)
    # No accuracy, significance, slope, bias, or outcome gate is applied.
    # PM and MINOBS only break an exact correlation tie.
    ranked <- eligible[order(-correlation[eligible], pm[eligible],
                             min_obs[eligible], as.character(group$pm_label[eligible]))]
    best <- group[ranked[[1]], , drop = FALSE]
    fit_path <- file.path(
      cfg$data_root, "data", "processed", "ALFA_K", ct,
      "ALFAK_fitnessLandscape", as.character(best$pm_label[[1]]),
      paste0("MINOBS_", as.integer(best$min_obs[[1]])), paste0(pid, ".Rds")
    )
    if (!file.exists(fit_path) || !is.finite(file.info(fit_path)$size) ||
        file.info(fit_path)$size <= 0L) {
      stop("Selected ALFA-K fit is missing: ", fit_path)
    }
    best
  })
  best <- do.call(rbind, picks)
  rownames(best) <- NULL
  selected_by_ct[[ct]] <- list(all_results = all_results, best_results = best)

  for (i in seq_len(nrow(best))) {
    row <- best[i, , drop = FALSE]
    pid <- as.character(row$PatientID[[1]])
    old <- original[original$cancer_type == ct & original$PatientID == pid, , drop = FALSE]
    sig <- as.logical(all_results$sig_flag[all_results$PatientID == pid])
    if (anyNA(sig)) stop("Invalid original sig_flag for ", ct, "/", pid)
    selected_rows[[length(selected_rows) + 1L]] <- data.frame(
      cancer_type = ct, PatientID = pid,
      pm_label = as.character(row$pm_label[[1]]),
      min_obs = as.integer(row$min_obs[[1]]),
      stringsAsFactors = FALSE
    )
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      cancer_type = ct, PatientID = pid,
      original_significant_candidates = sum(sig),
      original_selected = nrow(old) == 1L,
      original_pm_label = if (nrow(old)) as.character(old$pm_label[[1]]) else "",
      original_min_obs = if (nrow(old)) as.integer(old$min_obs[[1]]) else NA_integer_,
      correlation_pm_label = as.character(row$pm_label[[1]]),
      correlation_min_obs = as.integer(row$min_obs[[1]]),
      pearson_r = as.numeric(row$correlation[[1]]),
      cv_points = as.integer(row$total_points[[1]]),
      original_rule_passed = as.logical(row$sig_flag[[1]]),
      pearson_p = as.numeric(row$p_value[[1]]),
      permutation_p = as.numeric(row$huber_perm_p[[1]]),
      huber_slope = as.numeric(row$huber_slope[[1]]),
      signed_r2_huber_weighted = as.numeric(row$signed_r2_Huber_weighted[[1]]),
      stringsAsFactors = FALSE
    )
  }
}
selection <- do.call(rbind, selected_rows)
comparison <- do.call(rbind, comparison_rows)
if (nrow(selection) != nrow(expected) ||
    !setequal(paste(selection$cancer_type, selection$PatientID, sep = "\r"),
              expected_key) ||
    nrow(comparison) != nrow(expected)) {
  stop("Correlation-only selection does not cover all expected combinations")
}
if (sum(comparison$original_significant_candidates) != 10L ||
    sum(comparison$original_selected) != 5L ||
    sum(comparison$original_rule_passed) != 1L ||
    total_evaluated != 8820L) {
  stop("KFL comparison changed; review SELECTION_METHOD.md")
}

cohort_dir <- file.path(run_dir, "analysis", "results", "OptimizedParameters",
                        variant_cfg$workflow$analysis_id)
empty_abm <- data.frame(ct = character(), PatientID = character(),
                        pm_label = character(), min_obs = integer(),
                        stringsAsFactors = FALSE)
final <- list()
for (ct in names(selected_by_ct)) {
  best <- selected_by_ct[[ct]]$best_results
  tuple <- data.frame(
    ct = ct, PatientID = as.character(best$PatientID),
    pm_label = as.character(best$pm_label),
    min_obs = as.integer(best$min_obs),
    stringsAsFactors = FALSE
  )
  final[[ct]] <- list(kfl_best = tuple, ABM_best = empty_abm,
                      best_results = tuple)
  ct_cohort <- file.path(cohort_dir, ct)
  ct_root <- file.path(run_dir, "analysis", "results", "OptimizedParameters", ct)
  atomic_copy(file.path(source_results, analysis_id, ct, "correlation_results.csv"),
              file.path(ct_cohort, "correlation_results.csv"))
  atomic_copy(file.path(source_results, ct, "correlation_results.csv"),
              file.path(ct_root, "correlation_results.csv"))
  atomic_csv(best, file.path(ct_cohort, "best_correlation_results.csv"))
  atomic_csv(tuple, file.path(ct_cohort, "final_parameters_opt.csv"))
  atomic_csv(best, file.path(ct_root, "best_correlation_results.csv"))
  atomic_rds(selected_by_ct[[ct]], file.path(ct_root, "KFL_metrics.Rds"))
}
atomic_rds(selected_by_ct, file.path(cohort_dir, "parameters_opt.Rds"))
atomic_rds(final, file.path(cohort_dir, "final_parameters_opt.Rds"))
summary <- merge(expected, comparison,
                 by = c("cancer_type", "PatientID"), sort = FALSE)
summary <- summary[match(expected_key,
                         paste(summary$cancer_type, summary$PatientID, sep = "\r")), ,
                   drop = FALSE]
summary$time_start_day <- 0L
summary$time_end_day <- as.integer(unlist(cfg$patients[summary$PatientID]))
summary$time_unit <- "day"
summary$mapping_role <- ifelse(summary$high_cn == 6L, "primary", "sensitivity")
atomic_tsv(comparison, file.path(run_dir, "selection_comparison.tsv"))
atomic_tsv(summary, file.path(run_dir, "selection_summary.tsv"))
atomic_copy(file.path(cfg$project_root, "SELECTION_METHOD.md"),
            file.path(run_dir, "SELECTION_METHOD.md"))
provenance <- data.frame(
  key = c("method", "source_run_dir", "source_original_selection_sha256",
          "source_fit_audit_sha256", "source_analysis_id", "variant_analysis_id",
          "selected_combinations", "original_selected_combinations",
          "cart_git_sha", "pan_git_sha", "sif_sha256", "alfakr_commit",
          paste0("source_kfl_metrics_sha256_", names(selected_by_ct))),
  value = c("max_finite_pearson_correlation", source_dir,
            sha256(original_path), sha256(file.path(source_dir, "fit_audit.tsv")),
            analysis_id, variant_cfg$workflow$analysis_id,
            as.character(nrow(selection)), as.character(nrow(original)),
            trimws(system2("git", c("-C", cfg$project_root, "rev-parse", "HEAD"),
                           stdout = TRUE)),
            cfg$pan_commit, cfg$sif_sha256, cfg$alfakr_commit,
            vapply(names(selected_by_ct), function(ct) {
              sha256(file.path(source_results, ct, "KFL_metrics.Rds"))
            }, character(1))),
  stringsAsFactors = FALSE
)
atomic_tsv(provenance, file.path(run_dir, "selection_provenance.tsv"))

# Write the downstream controller's manifest last, after all parameter
# interfaces and provenance files have been created.
check <- readRDS(file.path(cohort_dir, "final_parameters_opt.Rds"))
if (!setequal(names(check), unique(expected$cancer_type)) ||
    sum(vapply(check, function(x) nrow(x$best_results), integer(1))) !=
      nrow(expected)) {
  stop("Final parameter RDS failed readback")
}
atomic_tsv(selection, file.path(run_dir, "manifests", "selected_samples.tsv"))
cat("Correlation-only selection complete: ", nrow(selection),
    " patient/mapping combinations; original-rule selections: ",
    nrow(original), ".\n", sep = "")
