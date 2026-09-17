args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: msr_reduce_first_minimum.R DOWNSTREAM_CONFIG VARIANT_RUN_DIR")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project <- dirname(dirname(normalizePath(script)))
source(file.path(project, "scripts", "common.R"))

config <- yaml::read_yaml(args[[1]])
run_dir <- normalizePath(args[[2]], mustWork = TRUE)
ct <- "CarT_high_cn_8"
if (!identical(basename(run_dir), "correlation_only") ||
    !identical(config$workflow$cart_selection_method,
               "max_finite_pearson_correlation")) {
  stop("This edge-case handler is limited to the Pearson-only run")
}
if (!identical(normalizePath(config$Data_path, mustWork = TRUE),
               file.path(run_dir, "analysis"))) {
  stop("Downstream Data_path does not match the variant run")
}
manifest <- utils::read.delim(
  file.path(run_dir, "manifests", "msr_samples.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
patient_info <- manifest[
  as.character(manifest$cancer_type) == ct,
  c("PatientID", "pm_label", "min_obs", "BKP"),
  drop = FALSE
]
if (nrow(patient_info) != 4L ||
    !setequal(as.character(patient_info$PatientID), c("P1", "P4", "P7", "P13"))) {
  stop("Unexpected high_cn_8 MSR patient manifest")
}

for (file in c(
  "code/packages/kflTopology/Ultis_MSR.R",
  "code/packages/kflVis/Vis_MSR.R",
  "code/packages/kflTopology/run_MSR.R"
)) {
  source(file.path(config$Repo_path, file), local = .GlobalEnv)
}
cache_path <- file.path(config$Data_path, "results", "MSR", ct,
                        "steadyStatePredictions.Rds")
cache <- readRDS(cache_path)
if (!setequal(names(cache), as.character(patient_info$PatientID))) {
  stop("Steady-state cache patient set changed")
}
minimum_index <- vapply(cache, function(x) which.min(as.numeric(x[1L, ])),
                        integer(1))
if (minimum_index[["P7"]] != 1L ||
    any(minimum_index[names(minimum_index) != "P7"] < 2L)) {
  stop("The validated first-minimum edge case changed")
}

original_compute_ETs <- compute_ETs
compute_ETs <- function(p_seq, H_curve) {
  if (length(p_seq) == 1L && length(H_curve) == 1L &&
      is.finite(p_seq[[1]]) && is.finite(H_curve[[1]])) {
    # There is no pre-minimum interval from which to estimate a threshold.
    return(c(ET_reg = NA_real_, ET_elbow = NA_real_, ET_cp = NA_real_,
             ET_curve = NA_real_, ET_slope = NA_real_, ET_model = NA_real_))
  }
  original_compute_ETs(p_seq, H_curve)
}

result <- run_MSR(
  Cancer_type = ct,
  patient_info = patient_info,
  data_dir = file.path(config$Data_path, "data", "processed", "ALFA_K"),
  results_dir = file.path(config$Data_path, "results", "MSR")
)
output_dir <- file.path(config$Data_path, "results", "MSR", ct)
atomic_plot_pdf <- function(plot, path, width = 5, height = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  ok <- FALSE
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (!ok) unlink(tmp)
  }, add = TRUE)
  grDevices::pdf(tmp, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  if (!file.rename(tmp, path)) stop("Cannot install plot: ", path)
  ok <- TRUE
}
atomic_plot_pdf(plot_curve(result$entropy_curve),
                file.path(output_dir, "Entropy_curve.pdf"))
atomic_plot_pdf(plot_curve(result$ss_curve),
                file.path(output_dir, "SS_curve.pdf"))

thresholds <- utils::read.csv(file.path(output_dir, "ET_summary.csv"),
                              stringsAsFactors = FALSE)
p7 <- thresholds[thresholds$sample == "P7", , drop = FALSE]
et_columns <- c("ET_reg", "ET_elbow", "ET_cp", "ET_curve", "ET_slope", "ET_model")
if (nrow(p7) != 1L || !all(et_columns %in% names(p7)) ||
    !all(is.na(p7[1L, et_columns]))) {
  stop("P7 first-minimum thresholds were not recorded as undefined")
}
diagnostic <- data.frame(
  cancer_type = ct, PatientID = names(minimum_index),
  steady_state_minimum_index = as.integer(minimum_index),
  threshold_handling = ifelse(
    names(minimum_index) == "P7", "six_undefined_thresholds_recorded_NA",
    "unaltered_PANcanKFLs_computation"
  ),
  stringsAsFactors = FALSE
)
atomic_tsv(diagnostic, file.path(run_dir, "msr_edge_case.tsv"))
if (!file.copy(file.path(project, "MSR_EDGE_CASE.md"),
               file.path(run_dir, "MSR_EDGE_CASE.md"), overwrite = TRUE)) {
  stop("Cannot copy MSR_EDGE_CASE.md into the run directory")
}

status <- data.frame(
  stage = "msr_reduce", cancer_type = ct, PatientID = "",
  state = "COMPLETED",
  message = paste0("samples=", nrow(patient_info),
                   ";P7_first_minimum_ETs=NA;output_dir=", output_dir),
  job_id = Sys.getenv("SLURM_JOB_ID"), array_task_id = "",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
atomic_tsv(status, file.path(run_dir, "status",
                             paste0("msr_reduce__", ct, ".tsv")))
cat("MSR high_cn_8 reduction complete; P7 first-minimum thresholds are NA.\n")
