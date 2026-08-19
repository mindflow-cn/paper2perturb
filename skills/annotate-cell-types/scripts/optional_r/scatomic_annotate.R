#!/usr/bin/env Rscript
#
# scATOMIC pan-cancer/TME expression annotation (optional adapter).
#
# Usage:
#   Rscript scatomic_annotate.R input.h5ad output.h5ad
#
# Only invoked when --mode deep --enable-r-adapters are both set,
# and strategy_category is tumor_or_tme.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scatomic_annotate.R input.h5ad output.h5ad")
}

input_path  <- args[1]
output_path <- args[2]

cat(sprintf("scATOMIC annotation: input=%s, output=%s\n", input_path, output_path))

required_packages <- c("scATOMIC", "reticulate", "zellkonverter")
missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing) > 0) {
  msg <- sprintf(
    "scATOMIC adapter skipped: missing R packages: %s. Install with:\n  devtools::install_github('abelson-lab/scATOMIC') for scATOMIC",
    paste(missing, collapse = ", ")
  )
  cat(msg, "\n", file = stderr())
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(scATOMIC)
})

# Load AnnData
cat("Loading AnnData...\n")
sce <- tryCatch({
  zellkonverter::readH5AD(input_path)
}, error = function(e) {
  cat(sprintf("Failed to read AnnData: %s\n", e$message), file = stderr())
  quit(status = 0)
})

# scATOMIC requires counts in the 'counts' assay
if (!"counts" %in% names(SummarizedExperiment::assays(sce))) {
  cat("Counts assay not found — cannot run scATOMIC\n", file = stderr())
  quit(status = 0)
}

cat("Running scATOMIC...\n")
pred <- tryCatch({
  scATOMIC::run_scATOMIC(sce)
}, error = function(e) {
  cat(sprintf("scATOMIC prediction failed: %s\n", e$message), file = stderr())
  quit(status = 0)
})

cat("Writing output...\n")
zellkonverter::writeH5AD(pred, output_path)

summary <- list(
  method = "scatomic",
  cells_labeled = ncol(pred)
)
if ("jsonlite" %in% installed.packages()[,"Package"]) {
  writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE), stdout())
}

cat("scATOMIC annotation complete.\n")
