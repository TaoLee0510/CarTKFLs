args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2L, 3L)) ||
    (length(args) == 3L && args[[3]] != "--dry-run")) {
  stop("Usage: select_r2_only.R CONFIG RUN_DIR [--dry-run]")
}
dry_run <- length(args) == 3L
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- normalizePath(args[[2]], mustWork = TRUE)
metric <- as.character(cfg$selection$r2_metric)
if (length(metric) != 1L || !metric %in% c("r2_unweighted", "r_squared", "R2R")) {
  stop("selection$r2_metric must be r2_unweighted, r_squared, or R2R")
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

if (identical(metric, "R2R")) {
  score <- vapply(seq_len(nrow(metrics)), function(i) {
    row <- metrics[i, , drop = FALSE]
    path <- file.path(fit_directory(cfg, row$pm_value, row$min_obs,
                                    row$PatientID), "xval.Rds")
    value <- readRDS(path)$R2R
    if (!is.numeric(value) || length(value) != 1L) {
      stop("Invalid ALFA-K R2R: ", path)
    }
    as.numeric(value)
  }, numeric(1))
} else {
  if (!metric %in% names(metrics)) stop("Missing R2 metric: ", metric)
  score <- suppressWarnings(as.numeric(metrics[[metric]]))
}
metrics$selection_r2 <- score
metrics$selection_r2_metric <- metric

best <- do.call(rbind, lapply(patients, function(patient) {
  rows <- metrics[metrics$PatientID == patient & is.finite(metrics$selection_r2),
                  , drop = FALSE]
  if (!nrow(rows)) stop("No finite ", metric, " candidate for ", patient)
  # Only R2 determines the ranking. The remaining keys resolve exact ties.
  order_idx <- order(-rows$selection_r2, rows$pm_value, -rows$min_obs)
  rows[order_idx[1L], , drop = FALSE]
}))
row.names(best) <- NULL
if (!identical(as.character(best$PatientID), patients)) {
  stop("R2 selection did not retain P16, P32, and P33")
}
if (dry_run) {
  print(best[, c("PatientID", "pm_label", "min_obs", "selection_r2",
                 "total_points"), drop = FALSE], row.names = FALSE)
  quit(save = "no", status = 0L)
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
source(file.path(cfg$pan_repo, "code/packages/kflInfer/parameters_optimize.R"))
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
  finite_r2_candidates = vapply(patients, function(patient) {
    sum(metrics$PatientID == patient & is.finite(metrics$selection_r2))
  }, integer(1)),
  significant_candidates_diagnostic = vapply(patients, function(patient) {
    sum(metrics$sig_flag[metrics$PatientID == patient] %in% TRUE)
  }, integer(1)),
  selection_status = "SELECTED_R2_ONLY",
  r2_metric = metric,
  best_r2 = best$selection_r2,
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
atomic_tsv(selected_rows, file.path(run_dir, "selection_summary.tsv"))
atomic_tsv(selected_rows[, c("cancer_type", "PatientID", "pm_label", "min_obs")],
           file.path(run_dir, "manifests", "selected_samples.tsv"))
revision <- data.frame(
  selected_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  selection_mode = "r2_only", r2_metric = metric,
  cartkfls_commit = trimws(system2(
    "git", c("-C", cfg$project_root, "rev-parse", "HEAD"), stdout = TRUE
  )),
  pancankfls_commit = cfg$pan_commit,
  sif_sha256 = cfg$sif_sha256,
  patient_count = nrow(selected_rows),
  downstream_state = "PENDING",
  stringsAsFactors = FALSE
)
atomic_tsv(revision, file.path(run_dir, "selection_revision.tsv"))
cat(sprintf("R2-only selection (%s): %d patients selected.\n", metric,
            nrow(selected_rows)))
