args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: summarize_selection.R CONFIG RUN_DIR")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
cfg <- read_config(args[[1]])
run_dir <- args[[2]]
selected_path <- file.path(run_dir, "manifests", "selected_samples.tsv")
selected <- utils::read.delim(selected_path, stringsAsFactors = FALSE, check.names = FALSE)
expected <- expand.grid(PatientID = names(cfg$patients),
                        high_cn = as.integer(unlist(cfg$high_cn)),
                        stringsAsFactors = FALSE)
expected$cancer_type <- vapply(expected$high_cn, cancer_type, character(1))
key <- paste(selected$cancer_type, selected$PatientID, sep = "\r")
expected_key <- paste(expected$cancer_type, expected$PatientID, sep = "\r")
if (nrow(selected) != nrow(expected) || anyDuplicated(key) ||
    !setequal(key, expected_key)) stop("KFL selection did not produce one row per patient/mapping")
selected$high_cn <- as.integer(sub("^CarT_high_cn_", "", selected$cancer_type))
selected$time_start_day <- 0L
selected$time_end_day <- vapply(selected$PatientID, function(x) as.integer(cfg$patients[[x]]), integer(1))
selected$time_unit <- "day"
selected$mapping_role <- ifelse(selected$high_cn == 6L, "primary", "sensitivity")
atomic_tsv(selected, file.path(run_dir, "selection_summary.tsv"))
cat("Selection summary passed for ", nrow(selected), " patient/mapping combinations.\n", sep = "")
