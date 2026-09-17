args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: diagnose_msr_thresholds.R STEADY_STATE_RDS")
predictions <- readRDS(args[[1]])
if (!is.list(predictions) || !length(predictions)) stop("Invalid steady-state cache")
for (patient in names(predictions)) {
  ss <- predictions[[patient]]
  if (!is.matrix(ss) || !nrow(ss) || !ncol(ss)) stop("Invalid matrix: ", patient)
  p <- suppressWarnings(as.numeric(colnames(ss)))
  curve <- as.numeric(ss[1L, ])
  idx <- which.min(curve)
  sub_p <- if (length(idx)) p[seq_len(idx)] else numeric()
  sub_curve <- if (length(idx)) curve[seq_len(idx)] else numeric()
  fit <- tryCatch(stats::lm(sub_curve ~ sub_p), error = function(e) e)
  slope <- if (inherits(fit, "error")) NA_real_ else unname(stats::coef(fit)[2L])
  cat(paste(
    patient, nrow(ss), ncol(ss), sum(is.finite(p)), sum(is.finite(curve)),
    if (length(idx)) idx else "NA", length(sub_curve),
    if (length(curve)) format(curve[1L], digits = 6) else "NA",
    if (length(curve)) format(tail(curve, 1L), digits = 6) else "NA",
    format(slope, digits = 6),
    sep = "\t"
  ), "\n", sep = "")
}
