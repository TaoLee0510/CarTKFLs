args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: audit_fit.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
tasks <- utils::read.delim(file.path(run_dir, "tasks_all.tsv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  outdir <- fit_directory(cfg, task$high_cn, task$pm, task$min_obs, task$patient)
  status_path <- file.path(run_dir, "status", sprintf("fit_%05d.tsv", task$task_id))
  status <- if (file.exists(status_path)) tryCatch(
    utils::read.delim(status_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  ) else NULL
  recorded <- if (!is.null(status) && nrow(status) == 1L) as.character(status$state) else "MISSING_STATUS"
  raw_valid <- validate_raw(outdir)
  flat_path <- flat_fit_path(outdir, task$patient)
  flat_valid <- if (file.exists(flat_path)) tryCatch({
    fit <- readRDS(flat_path)
    is.list(fit) && is.numeric(fit$fit_boot) && length(fit$fit_boot) > 0L &&
      is.data.frame(fit$xv_res) && all(c("f_est", "f_xv") %in% names(fit$xv_res))
  }, error = function(e) FALSE) else FALSE
  accepted <- (recorded == "COMPLETE" && raw_valid && flat_valid) ||
    (recorded %in% c("NO_CV", "MODEL_ERROR") && !flat_valid)
  data.frame(task, recorded_state = recorded, raw_valid = raw_valid,
             flat_valid = flat_valid, accepted = accepted,
             status_path = status_path, output_dir = outdir,
             stringsAsFactors = FALSE)
})
audit <- do.call(rbind, rows)
atomic_tsv(audit, file.path(run_dir, "fit_audit.tsv"))
if (any(!audit$accepted)) {
  stop(sum(!audit$accepted), " fit tasks have missing or inconsistent status; see fit_audit.tsv")
}

sample_rows <- list()
for (high_cn in as.integer(unlist(cfg$high_cn))) {
  for (patient in names(cfg$patients)) {
    subset <- audit[audit$high_cn == high_cn & audit$patient == patient, , drop = FALSE]
    valid_n <- sum(subset$recorded_state == "COMPLETE")
    sample_rows[[length(sample_rows) + 1L]] <- data.frame(
      cancer_type = cancer_type(high_cn), PatientID = patient,
      complete_raw = valid_n,
      model_error = sum(subset$recorded_state == "MODEL_ERROR"),
      no_cross_validation = sum(subset$recorded_state == "NO_CV"),
      stringsAsFactors = FALSE
    )
  }
}
summary <- do.call(rbind, sample_rows)
atomic_tsv(summary, file.path(run_dir, "fit_sample_summary.tsv"))
if (any(summary$complete_raw == 0L)) {
  stop("At least one patient/mapping has no evaluable fit; see fit_sample_summary.tsv")
}
atomic_tsv(summary[, c("cancer_type", "PatientID", "complete_raw")],
           file.path(run_dir, "manifests", "kfl_samples.tsv"))
cat(sprintf("Fit audit passed: %d complete; %d model errors; %d no-CV; %d samples.\n",
            sum(audit$recorded_state == "COMPLETE"),
            sum(audit$recorded_state == "MODEL_ERROR"),
            sum(audit$recorded_state == "NO_CV"), nrow(summary)))
