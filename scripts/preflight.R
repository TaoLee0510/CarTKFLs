args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: preflight.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]

if (!dir.exists(cfg$handoff_root) || !dir.exists(cfg$pan_repo) ||
    !file.exists(cfg$sif_path)) stop("Handoff, PANcanKFLs repo, or SIF is missing")
pan_sha <- trimws(system2("git", c("-C", cfg$pan_repo, "rev-parse", "HEAD"), stdout = TRUE))
if (!identical(pan_sha, cfg$pan_commit)) stop("PANcanKFLs commit differs: ", pan_sha)
sif_sha <- sha256(cfg$sif_path)
if (!identical(sif_sha, cfg$sif_sha256)) stop("SIF checksum differs: ", sif_sha)
handoff_audit <- utils::read.delim(file.path(cfg$handoff_root, "rds_consumer_validation.tsv"),
                                   stringsAsFactors = FALSE, check.names = FALSE)
for (field in c("patient", "high_cn", "schema_and_timing", "repository_consumer_contract",
                "prior_counts_and_states")) {
  if (!field %in% names(handoff_audit)) stop("Handoff audit missing ", field)
}

input_rows <- list()
for (high_cn in as.integer(unlist(cfg$high_cn))) {
  for (patient in names(cfg$patients)) {
    src <- source_input(cfg, high_cn, patient)
    dest <- consumer_input(cfg, high_cn, patient)
    expected_end <- as.integer(cfg$patients[[patient]])
    verified <- validate_input(src, patient, expected_end)
    expected <- cfg$expected_input[[patient]]
    if (!identical(as.integer(c(verified$n_karyotypes, verified$n_pre, verified$n_post)),
                   as.integer(c(expected$karyotypes, expected$pre_cells, expected$post_cells)))) {
      stop("Karyotype or retained-cell counts differ from handoff for ", patient, "/", high_cn)
    }
    check <- handoff_audit[handoff_audit$patient == patient &
                             handoff_audit$high_cn == high_cn, , drop = FALSE]
    if (nrow(check) != 1L || any(as.character(check[1, c(
      "schema_and_timing", "repository_consumer_contract", "prior_counts_and_states"
    )]) != "PASS")) stop("Handoff consumer validation failed for ", patient, "/", high_cn)
    src_sha <- sha256(src)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(dest)) {
      if (!file.copy(src, dest, overwrite = FALSE)) stop("Cannot stage ", src)
    }
    if (!identical(sha256(dest), src_sha)) stop("Staged input hash differs: ", dest)
    input_rows[[length(input_rows) + 1L]] <- data.frame(
      patient = patient, high_cn = high_cn, cancer_type = cancer_type(high_cn),
      time_start_day = 0L, time_end_day = expected_end,
      n_karyotypes = verified$n_karyotypes, n_pre = verified$n_pre,
      n_post = verified$n_post, sha256 = src_sha,
      source_path = src, staged_path = dest, stringsAsFactors = FALSE
    )
  }
}
atomic_tsv(do.call(rbind, input_rows), file.path(run_dir, "input_provenance.tsv"))

pm <- pm_values(cfg)
minobs <- as.integer(unlist(cfg$minobs))
if (length(pm) != 247L || !identical(minobs, c(5L, 10L, 20L))) {
  stop("PM/MINOBS grid differs from PANcanKFLs production grid")
}
rows <- list()
for (min_obs in minobs) {
  for (high_cn in as.integer(unlist(cfg$high_cn))) {
    for (patient in names(cfg$patients)) {
      rows[[length(rows) + 1L]] <- data.frame(
        patient = patient, high_cn = high_cn, cancer_type = cancer_type(high_cn),
        pm = pm, pm_label = vapply(pm, pm_label, character(1)),
        min_obs = min_obs, stringsAsFactors = FALSE
      )
    }
  }
}
tasks <- do.call(rbind, rows)
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[, c("task_id", "patient", "high_cn", "cancer_type", "pm", "pm_label", "min_obs")]
atomic_tsv(tasks, file.path(run_dir, "tasks_all.tsv"))
for (min_obs in minobs) {
  subset <- tasks[tasks$min_obs == min_obs, , drop = FALSE]
  atomic_tsv(subset, file.path(run_dir, paste0("tasks_MINOBS_", min_obs, ".tsv")))
}

pan_cfg_path <- file.path(cfg$pan_repo, "config", "config_hpc_hnsc_melanoma_ovarian_downstream.yaml")
downstream <- yaml::read_yaml(pan_cfg_path)
downstream$Repo_path <- cfg$pan_repo
downstream$Data_path <- cfg$data_root
downstream$workflow$analysis_id <- "GSE296419_CarT_days"
downstream$workflow$target_cancer_types <- vapply(as.integer(unlist(cfg$high_cn)), cancer_type, character(1))
downstream$workflow$expected_patient_ids <- stats::setNames(
  rep(list(names(cfg$patients)), length(cfg$high_cn)), downstream$workflow$target_cancer_types
)
downstream$workflow$excluded_samples <- list()
downstream$workflow$scientific_exclusions <- list()
downstream$workflow$alfak_minobs <- minobs
downstream$workflow$alfak_pm_grid <- list(explicit = as.numeric(cfg$pm_explicit),
  sequence = cfg$pm_sequence)
downstream$workflow$parameter_selection_mode <- "kfl_only"
downstream$workflow$run_expression_annotation <- FALSE
downstream$workflow$run_survival <- FALSE
downstream$workflow$run_abm <- FALSE
downstream$workflow$run_msr <- TRUE
downstream$alfak <- cfg$alfak
downstream$hpc$sif_path <- cfg$sif_path
downstream$datasets <- list()
downstream$MetaData <- list()
downstream$canonical_delivery <- NULL
downstream$save_task_rows <- FALSE
yaml::write_yaml(downstream, file.path(run_dir, "downstream.yaml"))

provenance <- data.frame(
  key = c("cart_git_sha", "pan_git_sha", "sif_sha256", "alfakr_commit",
          "handoff_root", "time_unit", "selection_mode", "grid_pm", "grid_minobs"),
  value = c(trimws(system2("git", c("-C", cfg$project_root, "rev-parse", "HEAD"), stdout = TRUE)),
            pan_sha, sif_sha, cfg$alfakr_commit, cfg$handoff_root,
            "day", "kfl_only_outcome_blind", length(pm), paste(minobs, collapse = ",")),
  stringsAsFactors = FALSE
)
atomic_tsv(provenance, file.path(run_dir, "provenance.tsv"))
cat("Preflight passed: ", nrow(tasks), " fit tasks; ", nrow(input_rows),
    " staged input RDS; SIF and PAN SHA verified.\n", sep = "")
