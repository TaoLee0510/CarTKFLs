script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "..", "common.R"))

base_read_config <- read_config
read_config <- function(path) {
  cfg <- base_read_config(path)
  if (is.null(cfg$cancer_type) || is.null(cfg$metadata_path) ||
      !identical(cfg$time_unit, "day")) stop("Incomplete myeloma configuration")
  if (!identical(sort(names(cfg$patients)), c("P16", "P32", "P33"))) {
    stop("Expected exactly P16, P32, and P33")
  }
  cfg
}

consumer_input <- function(cfg, patient) file.path(
  cfg$data_root, "data", "processed", "ALFA_K", cfg$cancer_type,
  "ALFA-K_inputs", paste0(patient, ".Rds")
)

fit_directory <- function(cfg, pm, minobs, patient) file.path(
  cfg$data_root, "data", "processed", "ALFA_K", cfg$cancer_type,
  "ALFAK_fitnessLandscape", pm_label(pm), paste0("MINOBS_", minobs), patient
)

numbat_matrix_path <- function(cfg, patient) {
  iteration <- as.integer(cfg$patients[[patient]]$numbat_iteration)
  file.path(cfg$handoff_root, patient, "Results",
            sprintf("GSE210079_%s_it%d_cnv_matrix.RData", patient, iteration))
}

expected_days <- function(cfg, patient) {
  as.numeric(vapply(cfg$patients[[patient]]$timepoints, function(x) x$day, numeric(1)))
}

validate_input <- function(path, patient, cfg) {
  yi <- readRDS(path)
  days <- expected_days(cfg, patient)
  if (!is.list(yi) || !is.data.frame(yi$x) || ncol(yi$x) != length(days) ||
      nrow(yi$x) < 2L || !identical(names(yi$x), as.character(days)) ||
      any(!is.finite(days)) || days[[1]] != 0 || is.unsorted(days, strictly = TRUE) ||
      !is.numeric(yi$dt) || length(yi$dt) != 1L || yi$dt != 1) {
    stop("ALFA-K input structure/timing mismatch for ", patient, ": ", path)
  }
  values <- as.matrix(yi$x)
  if (!is.numeric(values) || anyNA(values) || any(values < 0) ||
      any(values != round(values)) || any(colSums(values) <= 0)) {
    stop("Noninteger, negative, missing, or zero-depth counts: ", path)
  }
  states <- strsplit(rownames(values), ".", fixed = TRUE)
  if (length(states) != nrow(values) || any(lengths(states) != 22L) ||
      anyDuplicated(rownames(values)) ||
      any(!grepl("^[0-9]+(\\.[0-9]+){21}$", rownames(values)))) {
    stop("Karyotype key contract failed: ", path)
  }
  if (any(rownames(values) == paste(rep("2", 22), collapse = "."))) {
    stop("The excluded all-2 state is present: ", path)
  }
  list(yi = yi, n_karyotypes = nrow(values), depths = colSums(values))
}

fit_memory <- function(cfg, patient, minobs) {
  value <- cfg$hpc$fit_mem[[patient]][[paste0("MINOBS_", minobs)]]
  if (is.null(value) || !grepl("^[0-9]+G$", value)) stop("Missing fit memory tier")
  value
}
