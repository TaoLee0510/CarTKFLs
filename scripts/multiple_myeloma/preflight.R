args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: preflight.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- normalizePath(args[[2]], mustWork = TRUE)

if (!dir.exists(cfg$handoff_root) || !file.exists(cfg$metadata_path) ||
    !dir.exists(cfg$pan_repo) || !file.exists(cfg$sif_path)) {
  stop("Numbat, metadata, PANcanKFLs, or SIF path is missing")
}
pan_sha <- trimws(system2("git", c("-C", cfg$pan_repo, "rev-parse", "HEAD"), stdout = TRUE))
if (!identical(pan_sha, cfg$pan_commit)) stop("PANcanKFLs commit differs: ", pan_sha)
sif_sha <- sha256(cfg$sif_path)
if (!identical(sif_sha, cfg$sif_sha256)) stop("SIF checksum differs: ", sif_sha)
metadata <- utils::read.delim(cfg$metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("GSM", "sample_title", "subject_id") %in% names(metadata))) {
  stop("GSE210079 metadata lacks GSM, sample_title, or subject_id")
}

pm <- pm_values(cfg)
minobs <- as.integer(unlist(cfg$minobs))
if (length(pm) != 247L || !identical(minobs, c(5L, 10L, 20L))) {
  stop("PM/MINOBS grid differs from PANcanKFLs production grid")
}

input_rows <- list()
support_rows <- list()
task_rows <- list()
source_rows <- list()
for (patient in names(cfg$patients)) {
  source_path <- numbat_matrix_path(cfg, patient)
  if (!file.exists(source_path)) stop("Missing final Numbat matrix: ", source_path)
  env <- new.env(parent = emptyenv())
  objects <- load(source_path, envir = env)
  if (!all(c("cnv_matrix", "anno", "clone_post") %in% objects)) {
    stop("Numbat RData lacks cnv_matrix, anno, or clone_post: ", source_path)
  }
  cn <- env$cnv_matrix
  anno <- env$anno
  clones <- env$clone_post
  if (!all(c("cell", "GSM", "raw_barcode", "sample_label", "patient", "timepoint") %in%
           names(anno)) ||
      !all(c("cell", "clone_opt", "compartment_opt") %in% names(clones)) ||
      !"cell" %in% names(cn)) stop("Numbat object schema mismatch: ", patient)
  if (anyDuplicated(anno$cell) || anyDuplicated(clones$cell) ||
      anyDuplicated(cn$cell) || !setequal(anno$cell, clones$cell) ||
      !setequal(anno$cell, cn$cell)) {
    stop("Numbat cell identities mismatch or contain duplicates: ", patient)
  }
  clones <- clones[match(anno$cell, clones$cell), , drop = FALSE]
  cn <- cn[match(anno$cell, cn$cell), , drop = FALSE]
  if (!identical(as.character(anno$cell), as.character(clones$cell)) ||
      !identical(as.character(anno$cell), as.character(cn$cell)) ||
      any(as.character(anno$patient) != patient) ||
      any(as.character(anno$cell) !=
          paste0(anno$GSM, "__", anno$raw_barcode))) {
    stop("Numbat annotation joins or patient labels mismatch: ", patient)
  }
  if (anyNA(clones$compartment_opt) ||
      !all(clones$compartment_opt %in% c("tumor", "normal"))) {
    stop("Unexpected Numbat compartment: ", patient)
  }

  timepoints <- cfg$patients[[patient]]$timepoints
  labels <- names(timepoints)
  days <- expected_days(cfg, patient)
  if (length(labels) < 2L || days[[1]] != 0 ||
      is.unsorted(days, strictly = TRUE) ||
      !setequal(unique(as.character(anno$timepoint)), labels)) {
    stop("Timepoint configuration or Numbat annotation mismatch: ", patient)
  }
  for (label in labels) {
    spec <- timepoints[[label]]
    rows <- anno$timepoint == label
    if (!any(rows) || !identical(unique(as.character(anno$GSM[rows])),
                                 as.character(spec$gsm)) ||
        !identical(unique(as.character(anno$sample_label[rows])),
                   as.character(spec$sample_label))) {
      stop("Numbat timepoint/GSM/sample mismatch: ", patient, "/", label)
    }
    meta <- metadata[metadata$GSM == spec$gsm, , drop = FALSE]
    if (nrow(meta) != 1L ||
        !identical(as.character(meta$sample_title), as.character(spec$sample_label)) ||
        !identical(as.character(meta$subject_id),
                   paste0("Patient# ", sub("^P", "", patient)))) {
      stop("GSE210079 metadata mismatch: ", patient, "/", label)
    }
    observed <- sum(rows & clones$compartment_opt == "tumor")
    if (observed != as.integer(spec$tumor_cells)) {
      stop("Numbat tumor count differs from configured source: ",
           patient, "/", label, " observed=", observed)
    }
  }

  segments <- grep("^[0-9]+:[0-9]+-[0-9]+$", names(cn), value = TRUE)
  if (!length(segments)) stop("No autosomal CNV segments: ", patient)
  parsed <- strcapture("^([0-9]+):([0-9]+)-([0-9]+)$", segments,
                       proto = list(chr = integer(), start = numeric(), end = numeric()))
  if (any(parsed$end < parsed$start) ||
      any(!parsed$chr %in% seq_len(22L))) stop("Invalid chromosome segment: ", patient)
  parsed$length <- parsed$end - parsed$start + 1
  selected <- vapply(seq_len(22L), function(chr) {
    choices <- which(parsed$chr == chr)
    if (length(choices)) segments[choices[which.max(parsed$length[choices])]] else NA_character_
  }, character(1))
  states <- matrix(2L, nrow(anno), 22L)
  for (j in seq_len(22L)) {
    if (is.na(selected[[j]])) next
    value <- as.character(cn[[selected[[j]]]])
    valid <- is.na(value) | value %in% c("", "neu", "loh", "amp", "bamp", "del")
    if (!all(valid)) stop("Unknown CNV state on chromosome ", j, " for ", patient)
    states[, j] <- ifelse(value == "amp", 3L,
                         ifelse(value == "bamp", 4L,
                                ifelse(value == "del", 1L, 2L)))
    states[is.na(states[, j]), j] <- 2L
  }
  karyotype <- apply(states, 1L, paste, collapse = ".")
  day <- days[match(as.character(anno$timepoint), labels)]
  tumor <- clones$compartment_opt == "tumor"
  diploid <- karyotype == paste(rep("2", 22L), collapse = ".")
  retained <- tumor & !diploid
  ledger <- data.frame(
    patient = patient, cell = anno$cell, raw_barcode = anno$raw_barcode,
    GSM = anno$GSM, sample_label = anno$sample_label,
    timepoint = anno$timepoint, day = day,
    compartment_opt = clones$compartment_opt, clone_opt = clones$clone_opt,
    p_cnv = clones$p_cnv, karyotype = karyotype,
    retained_non_diploid = retained, stringsAsFactors = FALSE
  )
  atomic_tsv(ledger, file.path(run_dir, "input_ledgers", paste0(patient, ".tsv")))
  keys <- sort(unique(karyotype[retained]))
  if (length(keys) < 2L) stop("Too few non-diploid karyotypes: ", patient)
  x <- matrix(0L, nrow = length(keys), ncol = length(days),
              dimnames = list(keys, as.character(days)))
  for (j in seq_along(days)) {
    cells <- retained & day == days[[j]]
    x[, j] <- tabulate(match(karyotype[cells], keys), nbins = length(keys))
    if (sum(x[, j]) < 1L) stop("Zero non-diploid depth: ", patient, "/", labels[[j]])
    input_rows[[length(input_rows) + 1L]] <- data.frame(
      patient = patient, timepoint = labels[[j]], day = days[[j]],
      GSM = as.character(timepoints[[j]]$gsm), assigned_tumor = sum(tumor & day == days[[j]]),
      excluded_diploid = sum(tumor & diploid & day == days[[j]]),
      retained_non_diploid = sum(x[, j]), stringsAsFactors = FALSE
    )
  }
  input_path <- consumer_input(cfg, patient)
  atomic_rds(list(x = as.data.frame(x, optional = TRUE),
                  pop.fitness = NULL, dt = 1), input_path)
  verified <- validate_input(input_path, patient, cfg)
  if (!identical(as.integer(verified$depths), as.integer(colSums(x)))) {
    stop("Staged ALFA-K depth mismatch: ", patient)
  }
  source_rows[[length(source_rows) + 1L]] <- data.frame(
    patient = patient, source = source_path, source_sha256 = sha256(source_path),
    input = input_path, input_sha256 = sha256(input_path),
    n_karyotypes = nrow(x), time_days = paste(days, collapse = ","),
    selected_chr_segments = paste(ifelse(is.na(selected), "diploid", selected), collapse = ";"),
    stringsAsFactors = FALSE
  )
  counts <- rowSums(x)
  for (m in minobs) {
    frequent <- counts >= m
    n_frequent <- sum(frequent)
    pm_max <- if (n_frequent) {
      max_cn <- max(rowSums(matrix(
        as.integer(unlist(strsplit(rownames(x)[frequent], ".", fixed = TRUE))),
        ncol = 22L, byrow = TRUE
      )))
      1 - 0.5^(1 / max_cn)
    } else NA_real_
    viable_pm <- if (n_frequent) pm[pm < pm_max] else numeric()
    support_rows[[length(support_rows) + 1L]] <- data.frame(
      patient = patient, min_obs = m, n_karyotypes = nrow(x),
      n_frequent = n_frequent, pm_max_strict = pm_max,
      supported_pm = length(viable_pm),
      excluded_pm = length(pm) - length(viable_pm),
      stringsAsFactors = FALSE
    )
    for (value in viable_pm) {
      task_rows[[length(task_rows) + 1L]] <- data.frame(
        patient = patient, cancer_type = cfg$cancer_type,
        pm = value, pm_label = pm_label(value), min_obs = m,
        stringsAsFactors = FALSE
      )
    }
  }
}
atomic_tsv(do.call(rbind, input_rows), file.path(run_dir, "input_counts.tsv"))
atomic_tsv(do.call(rbind, source_rows), file.path(run_dir, "input_provenance.tsv"))
atomic_tsv(do.call(rbind, support_rows), file.path(run_dir, "input_support.tsv"))
if (!length(task_rows)) stop("No viable ALFA-K grid tasks")
tasks <- do.call(rbind, task_rows)
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[, c("task_id", "patient", "cancer_type", "pm", "pm_label", "min_obs")]
atomic_tsv(tasks, file.path(run_dir, "tasks_all.tsv"))
resource_rows <- list()
for (patient in names(cfg$patients)) for (m in minobs) {
  subset <- tasks[tasks$patient == patient & tasks$min_obs == m, , drop = FALSE]
  atomic_tsv(subset, file.path(run_dir, sprintf("tasks_%s_MINOBS_%d.tsv", patient, m)))
  resource_rows[[length(resource_rows) + 1L]] <- data.frame(
    patient = patient, min_obs = m, tasks = nrow(subset),
    cpus = 1L, mem = fit_memory(cfg, patient, m),
    qos = cfg$hpc$qos, time = cfg$hpc$time, stringsAsFactors = FALSE
  )
}
atomic_tsv(do.call(rbind, resource_rows), file.path(run_dir, "fit_resources.tsv"))

pan_cfg_path <- file.path(cfg$pan_repo, "config",
                          "config_hpc_hnsc_melanoma_ovarian_downstream.yaml")
downstream <- yaml::read_yaml(pan_cfg_path)
downstream$Repo_path <- cfg$pan_repo
downstream$Data_path <- cfg$data_root
downstream$workflow$analysis_id <- "GSE210079_MultipleMyeloma_days"
downstream$workflow$target_cancer_types <- cfg$cancer_type
downstream$workflow$expected_patient_ids <- setNames(list(names(cfg$patients)),
                                                      cfg$cancer_type)
downstream$workflow$excluded_samples <- list()
downstream$workflow$scientific_exclusions <- list()
downstream$workflow$alfak_minobs <- minobs
downstream$workflow$alfak_pm_grid <- list(explicit = as.numeric(cfg$pm_explicit),
                                         sequence = cfg$pm_sequence)
downstream$workflow$parameter_selection_mode <- "kfl_only"
downstream$workflow$kfl_selection <- list(mode = "correlation_only",
                                          metric = "pearson")
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
          "metadata_path", "metadata_sha256", "numbat_root", "time_unit",
          "month3_nominal_days", "selection_mode", "grid_pm", "grid_minobs"),
  value = c(trimws(system2("git", c("-C", cfg$project_root, "rev-parse", "HEAD"),
                            stdout = TRUE)), pan_sha, sif_sha, cfg$alfakr_commit,
            cfg$metadata_path, sha256(cfg$metadata_path), cfg$handoff_root, "day",
            cfg$month3_nominal_days,
            "correlation_only:pearson", length(pm),
            paste(minobs, collapse = ",")), stringsAsFactors = FALSE
)
atomic_tsv(provenance, file.path(run_dir, "provenance.tsv"))
cat(sprintf("Preflight passed: %d patients, %d viable fit tasks.\n",
            length(cfg$patients), nrow(tasks)))
