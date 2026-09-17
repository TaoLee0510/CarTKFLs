args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: summarize_selection.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
selected <- utils::read.delim(file.path(run_dir, "manifests", "selected_samples.tsv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cancer_type", "PatientID", "pm_label", "min_obs")
if (!all(required %in% names(selected))) stop("KFL selection schema mismatch")
if (anyDuplicated(selected$PatientID) || !all(selected$PatientID %in% names(cfg$patients)) ||
    any(selected$cancer_type != cfg$cancer_type)) {
  stop("KFL selection contains duplicate or unexpected patient")
}
fit_summary <- utils::read.delim(file.path(run_dir, "fit_sample_summary.tsv"),
                                 stringsAsFactors = FALSE)
analysis_id <- yaml::read_yaml(file.path(run_dir, "downstream.yaml"))$workflow$analysis_id
metrics_path <- file.path(run_dir, "analysis", "results", "OptimizedParameters",
                          analysis_id, cfg$cancer_type, "correlation_results.csv")
metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("PatientID", "sig_flag") %in% names(metrics))) {
  stop("KFL metrics lack PatientID or sig_flag")
}
expected <- data.frame(cancer_type = cfg$cancer_type,
                       PatientID = names(cfg$patients), stringsAsFactors = FALSE)
if (!setequal(unique(as.character(metrics$PatientID)),
              as.character(fit_summary$PatientID[fit_summary$complete_raw > 0]))) {
  stop("KFL metric coverage differs from evaluable fit manifest")
}
sig <- as.logical(metrics$sig_flag)
if (anyNA(sig)) stop("Invalid KFL sig_flag")
coverage <- expected
coverage$complete_fits <- fit_summary$complete_raw[match(expected$PatientID,
                                                          fit_summary$PatientID)]
coverage$significant_candidates <- vapply(expected$PatientID, function(patient) {
  sum(sig[metrics$PatientID == patient])
}, integer(1))
selected_index <- match(expected$PatientID, selected$PatientID)
if (any(coverage$significant_candidates[!is.na(selected_index)] < 1L) ||
    any(coverage$significant_candidates[is.na(selected_index)] > 0L)) {
  stop("Selection does not match significant KFL candidate coverage")
}
coverage$selection_status <- ifelse(
  !is.na(selected_index), "SELECTED",
  ifelse(coverage$complete_fits == 0L, "NO_EVALUABLE_FIT",
         "NO_SIGNIFICANT_KFL_CANDIDATE")
)
coverage$pm_label <- ifelse(is.na(selected_index), "", selected$pm_label[selected_index])
coverage$min_obs <- ifelse(is.na(selected_index), NA_integer_,
                           as.integer(selected$min_obs[selected_index]))
atomic_tsv(coverage, file.path(run_dir, "selection_coverage.tsv"))
atomic_tsv(coverage[is.na(selected_index), , drop = FALSE],
           file.path(run_dir, "selection_exclusions.tsv"))
selected$time_start_day <- 0L
selected$time_end_day <- vapply(selected$PatientID, function(patient) {
  tail(expected_days(cfg, patient), 1L)
}, numeric(1))
selected$time_unit <- "day"
atomic_tsv(selected, file.path(run_dir, "selection_summary.tsv"))
cat(sprintf("Selection: %d selected; %d excluded.\n",
            nrow(selected), sum(is.na(selected_index))))
