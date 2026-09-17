args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: fit_one.R CONFIG RUN_DIR TASK_FILE ROW_NUMBER")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
tasks <- utils::read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
row_num <- as.integer(args[[4]])
if (!is.finite(row_num) || row_num < 1L || row_num > nrow(tasks)) stop("Invalid task row")
task <- tasks[row_num, , drop = FALSE]
patient <- as.character(task$patient)
high_cn <- as.integer(task$high_cn)
pm <- as.numeric(task$pm)
minobs <- as.integer(task$min_obs)
task_id <- as.integer(task$task_id)
input <- consumer_input(cfg, high_cn, patient)
outdir <- fit_directory(cfg, high_cn, pm, minobs, patient)
status_path <- file.path(run_dir, "status", sprintf("fit_%05d.tsv", task_id))

write_status <- function(state, message = "") {
  message <- gsub("[\r\n\t]+", " ", as.character(message))
  atomic_tsv(data.frame(
    task_id = task_id, patient = patient, high_cn = high_cn,
    cancer_type = cancer_type(high_cn), pm = pm, pm_label = pm_label(pm),
    min_obs = minobs, state = state, message = message,
    job_id = Sys.getenv("SLURM_JOB_ID"),
    array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  ), status_path)
}

if (file.exists(status_path)) {
  previous <- utils::read.delim(status_path, stringsAsFactors = FALSE)
  if (nrow(previous) == 1L && previous$state == "COMPLETE" &&
      validate_raw(outdir) && file.exists(flat_fit_path(outdir, patient))) {
    cat("Existing complete task ", task_id, " retained\n", sep = "")
    quit(save = "no", status = 0L)
  }
  stop("Task already has a status record; no automatic retry: ", status_path)
}

verified <- validate_input(input, patient, cfg$patients[[patient]])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260917L + task_id)
state <- ""
failure <- ""
tryCatch({
  if (!validate_raw(outdir)) {
    alfakR::alfak(
      yi = verified$yi, outdir = outdir, passage_times = NULL,
      minobs = minobs, nboot = as.integer(cfg$alfak$nboot),
      n0 = as.numeric(cfg$alfak$n0), nb = as.numeric(cfg$alfak$nb),
      pm = pm, correct_efflux = isTRUE(cfg$alfak$correct_efflux),
      allow_noninteger_counts = FALSE
    )
  }
  if (!validate_raw(outdir)) stop("ALFA-K did not produce a valid raw quartet")
  extract_patient_fit(outdir, patient, minobs)
  state <- "COMPLETE"
}, error = function(e) {
  failure <<- conditionMessage(e)
  state <<- if (validate_raw(outdir) && grepl("cross-validation|xval", failure,
                                              ignore.case = TRUE)) "NO_CV" else "MODEL_ERROR"
})
if (state != "COMPLETE") {
  writeLines(c(
    paste0("task_id=", task_id), paste0("patient=", patient),
    paste0("high_cn=", high_cn), paste0("pm=", pm),
    paste0("minobs=", minobs), paste0("state=", state),
    paste0("message=", failure)
  ), file.path(outdir, "alfak_failed.txt"))
}
write_status(state, if (state == "COMPLETE") "" else failure)
cat(sprintf("task=%d patient=%s high_cn=%d pm=%s minobs=%d state=%s %s\n",
            task_id, patient, high_cn, pm_label(pm), minobs, state, failure))
if (state == "MODEL_ERROR") quit(save = "no", status = 1L)
