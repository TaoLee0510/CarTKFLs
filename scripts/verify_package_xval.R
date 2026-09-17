args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: verify_package_xval.R BOOTSTRAP_RDS")
bootstrap <- readRDS(args[[1]])
xval <- get("xval", envir = asNamespace("alfakR"))(bootstrap)
if (!is.list(xval) || !is.numeric(xval$R2R) || length(xval$R2R) != 1L ||
    !is.data.frame(xval$xval_data) || nrow(xval$xval_data) < 2L ||
    !all(c("k", "observation", "prediction") %in% names(xval$xval_data))) {
  stop("Installed alfakR returned an invalid cross-validation result")
}
cat(sprintf("alfakR=%s xval_rows=%d R2R=%s\n",
            as.character(utils::packageVersion("alfakR")),
            nrow(xval$xval_data), format(xval$R2R, digits = 8)))
