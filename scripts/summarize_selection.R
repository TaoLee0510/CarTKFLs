args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: summarize_selection.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
selected_path <- file.path(run_dir, "manifests", "selected_samples.tsv")
selected <- utils::read.delim(selected_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cancer_type", "PatientID", "pm_label", "min_obs")
if (!all(required %in% names(selected))) {
  stop("KFL selection is missing required columns: ",
       paste(setdiff(required, names(selected)), collapse = ", "))
}
expected <- expand.grid(PatientID = names(cfg$patients),
                        high_cn = as.integer(unlist(cfg$high_cn)),
                        stringsAsFactors = FALSE)
expected$cancer_type <- vapply(expected$high_cn, cancer_type, character(1))
key <- paste(selected$cancer_type, selected$PatientID, sep = "\r")
expected_key <- paste(expected$cancer_type, expected$PatientID, sep = "\r")
if (!nrow(selected) || anyDuplicated(key) || !all(key %in% expected_key)) {
  stop("KFL selection has no rows, duplicate rows, or unexpected patient/mapping pairs")
}

# PANcanKFLs selects only combinations with a significant KFL candidate. Keep
# every expected combination visible without assigning parameters to failures.
analysis_id <- yaml::read_yaml(file.path(run_dir, "downstream.yaml"))$workflow$analysis_id
if (is.null(analysis_id) || length(analysis_id) != 1L || !nzchar(analysis_id)) {
  stop("Missing downstream analysis_id")
}
counts <- integer(nrow(expected))
for (ct in unique(expected$cancer_type)) {
  metrics_path <- file.path(run_dir, "analysis", "results", "OptimizedParameters",
                            analysis_id, ct, "correlation_results.csv")
  if (!file.exists(metrics_path)) stop("Missing KFL selection metrics: ", metrics_path)
  metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("PatientID", "sig_flag") %in% names(metrics))) {
    stop("KFL selection metrics lack PatientID or sig_flag: ", metrics_path)
  }
  if (!setequal(unique(as.character(metrics$PatientID)),
                expected$PatientID[expected$cancer_type == ct])) {
    stop("KFL selection metrics have incomplete patient coverage: ", metrics_path)
  }
  sig <- as.logical(metrics$sig_flag)
  if (anyNA(sig)) stop("KFL selection metrics contain invalid sig_flag: ", metrics_path)
  for (i in which(expected$cancer_type == ct)) {
    counts[[i]] <- sum(sig[metrics$PatientID == expected$PatientID[[i]]])
  }
}
selected_index <- match(expected_key, key)
if (any(counts[!is.na(selected_index)] < 1L) ||
    any(counts[is.na(selected_index)] > 0L)) {
  stop("Selected rows do not match significant KFL candidate coverage")
}
coverage <- expected
coverage$significant_candidates <- counts
coverage$selection_status <- ifelse(is.na(selected_index),
                                    "NO_SIGNIFICANT_KFL_CANDIDATE", "SELECTED")
coverage$pm_label <- ifelse(is.na(selected_index), "", selected$pm_label[selected_index])
coverage$min_obs <- ifelse(is.na(selected_index), NA_integer_,
                           as.integer(selected$min_obs[selected_index]))
atomic_tsv(coverage, file.path(run_dir, "selection_coverage.tsv"))
atomic_tsv(coverage[is.na(selected_index), , drop = FALSE],
           file.path(run_dir, "selection_exclusions.tsv"))

selected$high_cn <- as.integer(sub("^CarT_high_cn_", "", selected$cancer_type))
selected$time_start_day <- 0L
selected$time_end_day <- vapply(selected$PatientID, function(x) as.integer(cfg$patients[[x]]), integer(1))
selected$time_unit <- "day"
selected$mapping_role <- ifelse(selected$high_cn == 6L, "primary", "sensitivity")
atomic_tsv(selected, file.path(run_dir, "selection_summary.tsv"))
cat("Selection summary passed for ", nrow(selected), " selected and ",
    sum(is.na(selected_index)), " excluded patient/mapping combinations.\n", sep = "")
