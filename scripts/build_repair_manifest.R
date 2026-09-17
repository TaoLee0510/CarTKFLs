args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: build_repair_manifest.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
tasks <- utils::read.delim(file.path(run_dir, "tasks_all.tsv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
eligible <- logical(nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  path <- file.path(run_dir, "status", sprintf("fit_%05d.tsv", task$task_id))
  if (!file.exists(path)) next
  status <- tryCatch(utils::read.delim(path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(status) || nrow(status) != 1L ||
      status$state != "MODEL_ERROR" ||
      !grepl("subscript out of bounds", status$message, fixed = TRUE)) next
  outdir <- fit_directory(cfg, task$high_cn, task$pm, task$min_obs, task$patient)
  raw3 <- file.path(outdir, c("bootstrap_res.Rds", "landscape.Rds",
                               "landscape_posterior_samples.Rds"))
  eligible[[i]] <- all(file.exists(raw3)) && all(file.info(raw3)$size > 0) &&
    !file.exists(file.path(outdir, "xval.Rds"))
}
repair <- tasks[eligible, , drop = FALSE]
atomic_tsv(repair, file.path(run_dir, "repair_tasks.tsv"))
cat("Eligible xval-only repairs: ", nrow(repair), "\n", sep = "")
