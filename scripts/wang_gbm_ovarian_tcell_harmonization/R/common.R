suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

log_message <- function(...) {
  cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))
  flush.console()
}

read_analysis_config <- function(path) {
  if (!file.exists(path)) stop("Configuration file not found: ", path)
  cfg <- yaml::read_yaml(path)
  required <- c("analysis_name", "inputs", "output_root", "sif_path", "required_hostname")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Configuration lacks: ", paste(missing, collapse = ", "))
  cfg
}

activate_task_library <- function(run_dir, cfg) {
  rel <- cfg$runtime$r_library_relpath %||% "runtime_library"
  lib <- file.path(run_dir, rel)
  if (!dir.exists(lib)) stop("Task R library not found: ", lib)
  .libPaths(unique(c(lib, .libPaths())))
  invisible(lib)
}

ensure_run_directories <- function(run_dir) {
  dirs <- c(
    "provenance", "audit", "annotated_objects", "cell_metadata",
    "sample_summaries", "qc_tables", "figures", "association_results",
    "report", "logs"
  )
  for (d in dirs) dir.create(file.path(run_dir, d), recursive = TRUE, showWarnings = FALSE)
  invisible(file.path(run_dir, dirs))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA"
  )
}

read_tsv <- function(path, ...) {
  utils::read.delim(
    path, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    stringsAsFactors = FALSE, ...
  )
}

safe_character <- function(x) {
  if (inherits(x, "factor")) x <- as.character(x)
  if (inherits(x, c("POSIXct", "POSIXlt", "Date"))) return(as.character(x))
  as.character(x)
}

short_values <- function(x, n = 12L) {
  x <- safe_character(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- unique(x)
  if (!length(x)) return("")
  paste(utils::head(x, n), collapse = " | ")
}

candidate_metadata_columns <- function(meta) {
  name_hit <- grepl(
    "cell|type|annot|cluster|sample|patient|orig|state|subtype|lineage|cnv|fitness|chrom|condition|group|treatment",
    colnames(meta), ignore.case = TRUE
  )
  n_unique <- vapply(meta, function(x) length(unique(x[!is.na(x)])), integer(1))
  which(name_hit | n_unique <= 100L)
}

metadata_profile <- function(meta) {
  data.frame(
    column = colnames(meta),
    class = vapply(meta, function(x) paste(class(x), collapse = ";"), character(1)),
    n_nonmissing = vapply(meta, function(x) sum(!is.na(x)), integer(1)),
    n_missing = vapply(meta, function(x) sum(is.na(x)), integer(1)),
    n_unique_nonmissing = vapply(meta, function(x) length(unique(x[!is.na(x)])), integer(1)),
    examples = vapply(meta, short_values, character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

metadata_value_counts <- function(meta, dataset) {
  idx <- candidate_metadata_columns(meta)
  out <- lapply(idx, function(j) {
    x <- safe_character(meta[[j]])
    x[is.na(x) | !nzchar(x)] <- "<NA_OR_EMPTY>"
    tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
    if (length(tab) > 250L) tab <- utils::head(tab, 250L)
    data.frame(
      dataset = dataset,
      column = colnames(meta)[j],
      value = names(tab),
      n_cells = as.integer(tab),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

object_profile <- function(obj, dataset, source_path) {
  assays <- tryCatch(Seurat::Assays(obj), error = function(e) character())
  reductions <- tryCatch(names(obj@reductions), error = function(e) character())
  data.frame(
    dataset = dataset,
    source_path = source_path,
    object_class = paste(class(obj), collapse = ";"),
    n_genes = nrow(obj),
    n_cells = ncol(obj),
    default_assay = tryCatch(Seurat::DefaultAssay(obj), error = function(e) NA_character_),
    assays = paste(assays, collapse = ";"),
    reductions = paste(reductions, collapse = ";"),
    n_metadata_columns = ncol(obj@meta.data),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
