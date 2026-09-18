args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2L, 3L)) ||
    (length(args) == 3L && args[[3]] != "--dry-run")) {
  stop("Usage: select_correlation_only.R CONFIG RUN_DIR [--dry-run]")
}
dry_run <- length(args) == 3L
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- normalizePath(args[[2]], mustWork = TRUE)
metric <- as.character(cfg$selection$metric)
if (!identical(metric, "correlation")) {
  stop("selection$metric must be correlation (Pearson)")
}
downstream <- yaml::read_yaml(file.path(run_dir, "downstream.yaml"))
analysis_id <- downstream$workflow$analysis_id
if (is.null(analysis_id) || !nzchar(analysis_id) ||
    !identical(downstream$Data_path, cfg$data_root)) {
  stop("Downstream configuration differs from this run")
}
manifest <- utils::read.delim(file.path(run_dir, "manifests", "kfl_samples.tsv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
patients <- names(cfg$patients)
if (!setequal(manifest$PatientID, patients) || anyDuplicated(manifest$PatientID) ||
    any(manifest$cancer_type != cfg$cancer_type)) {
  stop("KFL sample manifest does not cover exactly the three patients")
}

read_metrics <- function(patient) {
  path <- file.path(run_dir, "sample_results", "kfl", cfg$cancer_type,
                    paste0(patient, ".Rds"))
  entry <- readRDS(path)
  value <- entry$all_results
  needed <- c("PatientID", "pm_label", "pm_value", "min_obs", "total_points")
  if (!is.data.frame(value) || !nrow(value) ||
      !all(needed %in% names(value)) || any(value$PatientID != patient)) {
    stop("Invalid KFL metrics for ", patient, ": ", path)
  }
  value
}
metrics <- do.call(rbind, lapply(patients, read_metrics))
row.names(metrics) <- NULL
expected_per_patient <- length(pm_values(cfg)) * length(cfg$minobs)
if (nrow(metrics) != expected_per_patient * length(patients) ||
    any(table(metrics$PatientID) != expected_per_patient) ||
    any(!is.finite(metrics$pm_value)) ||
    anyNA(metrics$min_obs) ||
    any(metrics$pm_label != pm_label(metrics$pm_value)) ||
    anyDuplicated(metrics[, c("PatientID", "pm_label", "min_obs")])) {
  stop("KFL metric grid is incomplete or inconsistent")
}
expected_grid <- expand.grid(
  PatientID = patients, pm_label = pm_label(pm_values(cfg)),
  min_obs = as.integer(cfg$minobs), stringsAsFactors = FALSE
)
grid_key <- function(x) paste(x$PatientID, x$pm_label, x$min_obs, sep = "\t")
if (!setequal(grid_key(metrics), grid_key(expected_grid))) {
  stop("KFL metric grid differs from the configured combinations")
}

needed_metrics <- c("correlation", "p_value", "r2_unweighted", "sig_flag",
                    "huber_perm_p", "huber_slope", "intercept_rel",
                    "outlier_ratio", "bias_score", "sd_score", "post_score")
if (!all(needed_metrics %in% names(metrics))) {
  stop("KFL metrics lack Pearson or strict-selection diagnostic fields")
}
if (anyNA(metrics$sig_flag)) stop("KFL significance flags contain NA")
metrics$selection_pearson <- suppressWarnings(as.numeric(metrics$correlation))
metrics$selection_method <- "max_finite_pearson_correlation"

# Reconstruct the pinned PANcanKFLs selection in scratch space so the report
# compares against that exact implementation without changing prior outputs.
source(file.path(cfg$pan_repo, "code/packages/kflInfer/parameters_optimize.R"))
reconstruct_strict <- function() {
  scratch <- tempfile("myeloma_strict_selection_")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  strict_cfg <- downstream
  strict_cfg$Data_path <- scratch
  strict <- KFLAccuracyStabilityOpt(
    stats::setNames(list(list(all_results = metrics)), cfg$cancer_type),
    config = strict_cfg, out_subdir = "OptimizedParameters"
  )
  strict[[cfg$cancer_type]]$best_results
}
strict_best <- reconstruct_strict()
if (!is.data.frame(strict_best) || nrow(strict_best) != 1L ||
    strict_best$PatientID != "P32" ||
    strict_best$pm_label != "pm_0.0001" || strict_best$min_obs != 20L) {
  stop("Original PANcanKFLs selection changed; review the method report")
}

best <- do.call(rbind, lapply(patients, function(patient) {
  rows <- metrics[metrics$PatientID == patient & is.finite(metrics$selection_pearson),
                  , drop = FALSE]
  if (!nrow(rows)) stop("No finite Pearson correlation for ", patient)
  # Only Pearson r determines the ranking. The other keys resolve exact ties.
  order_idx <- order(-rows$selection_pearson, rows$pm_value, rows$min_obs,
                     rows$pm_label)
  winner <- rows[order_idx[1L], , drop = FALSE]
  fit_path <- flat_fit_path(fit_directory(
    cfg, winner$pm_value, winner$min_obs, patient
  ), patient)
  if (!file.exists(fit_path) || file.info(fit_path)$size <= 0L) {
    stop("Selected ALFA-K fit is missing: ", fit_path)
  }
  winner
}))
row.names(best) <- NULL
if (!identical(as.character(best$PatientID), patients)) {
  stop("Pearson selection did not retain P16, P32, and P33")
}
if (any(best$sig_flag %in% TRUE)) {
  stop("A Pearson winner now passes the original rule; review method report")
}
comparison <- do.call(rbind, lapply(patients, function(patient) {
  winner <- best[best$PatientID == patient, , drop = FALSE]
  original <- strict_best[strict_best$PatientID == patient, , drop = FALSE]
  data.frame(
    cancer_type = cfg$cancer_type, PatientID = patient,
    completed_fits = as.integer(manifest$complete_raw[match(patient, manifest$PatientID)]),
    original_significant_candidates = sum(metrics$sig_flag[metrics$PatientID == patient] %in% TRUE),
    original_selected = nrow(original) == 1L,
    original_pm_label = if (nrow(original)) as.character(original$pm_label) else "",
    original_min_obs = if (nrow(original)) as.integer(original$min_obs) else NA_integer_,
    pearson_pm_label = as.character(winner$pm_label),
    pearson_min_obs = as.integer(winner$min_obs),
    pearson_r = as.numeric(winner$selection_pearson),
    pearson_p = as.numeric(winner$p_value),
    cv_points = as.integer(winner$total_points),
    predictive_r2 = as.numeric(winner$r2_unweighted),
    original_rule_passed = as.logical(winner$sig_flag),
    permutation_p = as.numeric(winner$huber_perm_p),
    huber_slope = as.numeric(winner$huber_slope),
    intercept_rel = as.numeric(winner$intercept_rel),
    outlier_ratio = as.numeric(winner$outlier_ratio),
    bias_score = as.numeric(winner$bias_score),
    stringsAsFactors = FALSE
  )
}))
if (!identical(as.integer(comparison$original_significant_candidates),
               c(0L, 6L, 0L))) {
  stop("Original significant-candidate counts changed; review method report")
}
if (dry_run) {
  print(best[, c("PatientID", "pm_label", "min_obs", "selection_pearson",
                 "r2_unweighted", "total_points", "sig_flag"), drop = FALSE],
        row.names = FALSE)
  quit(save = "no", status = 0L)
}

downstream$workflow$kfl_selection <- list(
  mode = "correlation_only", metric = "pearson",
  ranking = "max_finite_correlation"
)
config_tmp <- tempfile(".downstream_", tmpdir = run_dir)
yaml::write_yaml(downstream, config_tmp)
if (!file.rename(config_tmp, file.path(run_dir, "downstream.yaml"))) {
  stop("Cannot install Pearson downstream configuration")
}

cohort_subdir <- file.path("OptimizedParameters", analysis_id)
cohort_dir <- file.path(cfg$data_root, "results", cohort_subdir)
ct_dir <- file.path(cohort_dir, cfg$cancer_type)
legacy_ct_dir <- file.path(cfg$data_root, "results", "OptimizedParameters",
                           cfg$cancer_type)
dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legacy_ct_dir, recursive = TRUE, showWarnings = FALSE)
atomic_csv <- function(value, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(value, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Cannot install ", path)
}
selected <- stats::setNames(list(list(all_results = metrics, best_results = best)),
                            cfg$cancer_type)
final <- parameters_optimize(config = downstream, kfl_results = selected,
                             out_subdir = cohort_subdir, selection_mode = "kfl_only")
selected_rows <- data.frame(cancer_type = cfg$cancer_type,
                            PatientID = best$PatientID,
                            pm_label = best$pm_label,
                            min_obs = as.integer(best$min_obs),
                            stringsAsFactors = FALSE)
if (!setequal(final[[cfg$cancer_type]]$best_results$PatientID, patients)) {
  stop("Final parameter structure does not contain all patients")
}
atomic_rds(selected[[cfg$cancer_type]], file.path(ct_dir, "KFL_metrics.Rds"))
atomic_rds(selected[[cfg$cancer_type]], file.path(legacy_ct_dir, "KFL_metrics.Rds"))
atomic_rds(selected, file.path(cohort_dir, "parameters_opt.Rds"))
for (out_dir in c(ct_dir, legacy_ct_dir)) {
  atomic_csv(metrics, file.path(out_dir, "correlation_results.csv"))
  atomic_csv(best, file.path(out_dir, "best_correlation_results.csv"))
}

coverage <- data.frame(
  cancer_type = cfg$cancer_type,
  PatientID = patients,
  complete_fits = as.integer(manifest$complete_raw[match(patients, manifest$PatientID)]),
  finite_correlation_candidates = vapply(patients, function(patient) {
    sum(metrics$PatientID == patient & is.finite(metrics$selection_pearson))
  }, integer(1)),
  significant_candidates_diagnostic = comparison$original_significant_candidates,
  selection_status = "SELECTED_PEARSON_ONLY",
  selection_metric = "pearson_correlation",
  best_correlation = best$selection_pearson,
  predictive_r2_diagnostic = best$r2_unweighted,
  total_points = as.integer(best$total_points),
  pm_label = best$pm_label,
  min_obs = as.integer(best$min_obs),
  stringsAsFactors = FALSE
)
selected_rows$time_start_day <- 0L
selected_rows$time_end_day <- vapply(selected_rows$PatientID, function(patient) {
  tail(expected_days(cfg, patient), 1L)
}, numeric(1))
selected_rows$time_unit <- "day"
atomic_tsv(coverage, file.path(run_dir, "selection_coverage.tsv"))
atomic_tsv(coverage[FALSE, , drop = FALSE], file.path(run_dir, "selection_exclusions.tsv"))
atomic_tsv(comparison, file.path(run_dir, "selection_comparison.tsv"))
atomic_tsv(selected_rows, file.path(run_dir, "selection_summary.tsv"))
atomic_tsv(selected_rows[, c("cancer_type", "PatientID", "pm_label", "min_obs")],
           file.path(run_dir, "manifests", "selected_samples.tsv"))
note_source <- file.path(cfg$project_root, "GSE210079_SELECTION_METHOD.md")
note_tmp <- tempfile(".SELECTION_METHOD_", tmpdir = run_dir)
if (!file.copy(note_source, note_tmp, overwrite = TRUE) ||
    !file.rename(note_tmp, file.path(run_dir, "SELECTION_METHOD.md"))) {
  stop("Cannot install the selection method note")
}
revision <- data.frame(
  selected_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  selection_mode = "correlation_only", selection_metric = "pearson",
  cartkfls_commit = trimws(system2(
    "git", c("-C", cfg$project_root, "rev-parse", "HEAD"), stdout = TRUE
  )),
  pancankfls_commit = cfg$pan_commit,
  sif_sha256 = cfg$sif_sha256,
  patient_count = nrow(selected_rows),
  downstream_state_at_selection = "PENDING",
  stringsAsFactors = FALSE
)
atomic_tsv(revision, file.path(run_dir, "selection_revision.tsv"))
cat(sprintf("Pearson-only selection: %d patients selected; strict rule selected %d.\n",
            nrow(selected_rows), nrow(strict_best)))
