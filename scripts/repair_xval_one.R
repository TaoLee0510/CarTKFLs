args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: repair_xval_one.R CONFIG RUN_DIR MANIFEST ROW_NUMBER")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
source(file.path(dirname(normalizePath(script)), "xval_compat.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
tasks <- utils::read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
index <- as.integer(args[[4]])
if (!is.finite(index) || index < 1L || index > nrow(tasks)) stop("Invalid repair row")
task <- tasks[index, , drop = FALSE]
outdir <- fit_directory(cfg, task$high_cn, task$pm, task$min_obs, task$patient)
status_path <- file.path(run_dir, "status", sprintf("fit_%05d.tsv", task$task_id))
status <- utils::read.delim(status_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(status) != 1L || status$state != "MODEL_ERROR" ||
    !grepl("subscript out of bounds", status$message, fixed = TRUE)) {
  stop("Repair target no longer has the audited dev2 error: ", status_path)
}
backup <- sub("[.]tsv$", ".before_xval_repair.tsv", status_path)
if (file.exists(backup)) stop("Repair was already attempted: ", backup)
raw3 <- file.path(outdir, c("bootstrap_res.Rds", "landscape.Rds",
                             "landscape_posterior_samples.Rds"))
if (!all(file.exists(raw3)) || any(file.info(raw3)$size <= 0) ||
    file.exists(file.path(outdir, "xval.Rds"))) {
  stop("Expected three complete raw files and no xval file: ", outdir)
}
set.seed(20260917L + as.integer(task$task_id))
xval <- xval_dev2_compat(readRDS(file.path(outdir, "bootstrap_res.Rds")))
atomic_rds(xval, file.path(outdir, "xval.Rds"))
extract_patient_fit(outdir, task$patient, task$min_obs)
if (!file.copy(status_path, backup, overwrite = FALSE)) stop("Cannot preserve prior status")
status$state <- "COMPLETE"
status$message <- "xval repaired from saved bootstrap; fitted outputs reused; see prior status backup"
status$repair_job_id <- Sys.getenv("SLURM_JOB_ID")
status$repair_timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
atomic_tsv(status, status_path)
writeLines(status$message, file.path(outdir, "xval_recovered.txt"))
cat(sprintf("Recovered task %d %s high_cn=%d %s MINOBS=%d R2R=%s\n",
            task$task_id, task$patient, task$high_cn, task$pm_label,
            task$min_obs, as.character(xval$R2R)))
