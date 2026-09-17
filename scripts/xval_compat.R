# The pinned alfakR dev2 SIF (source commit c5db5d8) names the second
# cross-validation column according to the Krig prediction object in some
# folds, but build_xval_result indexes it by the hard-coded name "est_f".
# Keep the installed fitting and cross-validation algorithm unchanged, except
# for indexing that second column by position. Fail closed if the pinned
# function body no longer contains exactly the audited expression.
xval_dev2_compat <- function(fq_boot) {
  original <- get("xval", envir = asNamespace("alfakR"))
  rendered <- deparse(body(original), width.cutoff = 500L)
  old <- 'tmp[, "est_f"]'
  replacement <- 'tmp[, 2L]'
  matches <- grepl(old, rendered, fixed = TRUE)
  if (sum(matches) != 1L) {
    stop("Installed alfakR xval body differs from the audited dev2 compatibility target")
  }
  rendered[matches] <- sub(old, replacement, rendered[matches], fixed = TRUE)
  patched <- original
  body(patched) <- parse(text = paste(rendered, collapse = "\n"))[[1L]]
  result <- patched(fq_boot)
  if (!is.list(result) || !is.numeric(result$R2R) || length(result$R2R) != 1L ||
      !is.data.frame(result$xval_data) ||
      !all(c("k", "observation", "prediction") %in% names(result$xval_data))) {
    stop("Compatibility xval result does not satisfy the dev2 output contract")
  }
  result
}
