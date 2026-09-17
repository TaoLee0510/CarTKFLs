`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is required")
  cfg <- yaml::read_yaml(normalizePath(path, mustWork = TRUE))
  for (field in c("project_root", "data_root", "handoff_root", "pan_repo", "sif_path")) {
    if (is.null(cfg[[field]]) || !nzchar(cfg[[field]])) stop("Missing config field: ", field)
  }
  cfg
}

pm_values <- function(cfg) {
  sort(unique(round(c(
    as.numeric(cfg$pm_explicit),
    seq(as.numeric(cfg$pm_sequence$start), as.numeric(cfg$pm_sequence$end),
        by = as.numeric(cfg$pm_sequence$by))
  ), 9)))
}

pm_label <- function(pm) {
  rendered <- format(as.numeric(pm), scientific = FALSE, trim = TRUE, digits = 15)
  rendered <- sub("0+$", "", rendered)
  rendered <- sub("[.]$", "", rendered)
  paste0("pm_", rendered)
}

cancer_type <- function(high_cn) paste0("CarT_high_cn_", as.integer(high_cn))

source_input <- function(cfg, high_cn, patient) file.path(
  cfg$handoff_root, "matrices", "majority_mean_exploratory",
  paste0("high_cn_", high_cn), patient, "ALFAK_input.Rds"
)

consumer_input <- function(cfg, high_cn, patient) file.path(
  cfg$data_root, "data", "processed", "ALFA_K", cancer_type(high_cn),
  "ALFA-K_inputs", paste0(patient, ".Rds")
)

fit_directory <- function(cfg, high_cn, pm, minobs, patient) file.path(
  cfg$data_root, "data", "processed", "ALFA_K", cancer_type(high_cn),
  "ALFAK_fitnessLandscape", pm_label(pm), paste0("MINOBS_", minobs), patient
)

flat_fit_path <- function(outdir, patient) file.path(dirname(outdir), paste0(patient, ".Rds"))

atomic_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(value, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Cannot install ", path)
  invisible(path)
}

atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("Cannot install ", path)
  invisible(path)
}

sha256 <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  output <- suppressWarnings(system2("sha256sum", path, stdout = TRUE, stderr = TRUE))
  if (length(output) != 1L || !grepl("^[0-9a-f]{64}  ", output)) {
    stop("sha256sum failed for ", path, ": ", paste(output, collapse = "; "))
  }
  substr(output, 1L, 64L)
}

validate_input <- function(path, patient, expected_end) {
  yi <- readRDS(path)
  if (!is.list(yi) || !is.data.frame(yi$x) || ncol(yi$x) != 2L ||
      nrow(yi$x) < 2L || !identical(names(yi$x), c("0", as.character(expected_end))) ||
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
    stop("The excluded all-2 fitting state is present in ", path)
  }
  list(yi = yi, n_karyotypes = nrow(values), n_pre = sum(values[, 1]),
       n_post = sum(values[, 2]))
}

required_raw <- c("bootstrap_res.Rds", "landscape.Rds",
                  "landscape_posterior_samples.Rds", "xval.Rds")

validate_raw <- function(outdir) {
  paths <- file.path(outdir, required_raw)
  if (!all(file.exists(paths)) || any(file.info(paths)$size <= 0)) return(FALSE)
  tryCatch({
    landscape <- readRDS(file.path(outdir, "landscape.Rds"))
    xval <- readRDS(file.path(outdir, "xval.Rds"))
    bootstrap <- readRDS(file.path(outdir, "bootstrap_res.Rds"))
    posterior <- readRDS(file.path(outdir, "landscape_posterior_samples.Rds"))
    is.data.frame(landscape) && nrow(landscape) > 0L &&
      all(c("k", "mean") %in% names(landscape)) &&
      is.list(xval) && !is.null(xval$R2R) &&
      !is.null(bootstrap) && !is.null(posterior)
  }, error = function(e) FALSE)
}

extract_patient_fit <- function(outdir, patient, minobs) {
  if (!validate_raw(outdir)) stop("Raw ALFA-K quartet invalid: ", outdir)
  landscape <- readRDS(file.path(outdir, "landscape.Rds"))
  xval <- readRDS(file.path(outdir, "xval.Rds"))
  fit_col <- if ("mean" %in% names(landscape)) "mean" else "median"
  if (!fit_col %in% names(landscape) || is.null(xval$tmp)) {
    stop("No evaluable cross-validation/fitness data: ", outdir)
  }
  f_mat <- as.data.frame(xval$tmp)
  if (ncol(f_mat) < 2L) stop("Cross-validation has fewer than two columns: ", outdir)
  names(f_mat)[1:2] <- c("f_est", "f_xv")
  fq <- if ("fq" %in% names(landscape)) {
    if (is.logical(landscape$fq)) landscape$fq else
      tolower(as.character(landscape$fq)) %in% c("true", "t", "1")
  } else rep(TRUE, nrow(landscape))
  xfq <- landscape[fq, , drop = FALSE]
  if (!nrow(xfq)) xfq <- landscape
  common <- intersect(rownames(f_mat), as.character(xfq$k))
  if (!length(common)) stop("No landscape/xval karyotype overlap: ", outdir)
  f_mat <- f_mat[common, c("f_est", "f_xv"), drop = FALSE]
  xfq <- xfq[match(common, xfq$k), , drop = FALSE]
  xfq$f_est <- suppressWarnings(as.numeric(f_mat$f_est))
  xfq$f_xv <- suppressWarnings(as.numeric(f_mat$f_xv))
  fit <- list(
    fit_boot = stats::setNames(as.numeric(landscape[[fit_col]]), as.character(landscape$k)),
    min_obs = as.integer(minobs), xv_res = xfq,
    vx_cor = tryCatch(stats::cor(xfq$f_est, xfq$f_xv, use = "complete"),
                      error = function(e) -Inf)
  )
  atomic_rds(fit, flat_fit_path(outdir, patient))
  invisible(fit)
}
