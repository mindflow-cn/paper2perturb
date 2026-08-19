#!/usr/bin/env Rscript
#
# SingleR reference-based cell-type annotation (optional adapter).
#
# Usage:
#   Rscript singler_annotate.R input.h5ad output.h5ad [--ref HumanPrimaryCellAtlasData]
#
# This adapter is only invoked when:
#   - --mode deep
#   - --enable-r-adapters
#   - required R packages are installed
#
# If dependencies are missing, exits with a clear message.
# The orchestrator falls back to Python methods on failure.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript singler_annotate.R input.h5ad output.h5ad [--ref reference_name]")
}

input_path  <- args[1]
output_path <- args[2]
ref_name    <- if (length(args) >= 4 && args[3] == "--ref") args[4] else "HumanPrimaryCellAtlasData"

cat(sprintf("SingleR annotation: input=%s, output=%s, ref=%s\n",
            input_path, output_path, ref_name))

# Check dependencies
required_packages <- c("SingleR", "celldex", "SingleCellExperiment", "zellkonverter", "reticulate")
missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing) > 0) {
  msg <- sprintf(
    "SingleR adapter skipped: missing R packages: %s. Install with:\n  BiocManager::install(c(%s))",
    paste(missing, collapse = ", "),
    paste(sprintf("'%s'", missing), collapse = ", ")
  )
  cat(msg, "\n", file = stderr())
  quit(status = 0)  # Exit 0 so orchestrator can fall back
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
})

# Load reference
cat(sprintf("Loading reference: %s\n", ref_name))
ref <- tryCatch({
  do.call(ref_name, list())
}, error = function(e) {
  cat(sprintf("Failed to load reference '%s': %s\n", ref_name, e$message),
      file = stderr())
  quit(status = 0)
})

# Load AnnData via zellkonverter
cat("Loading AnnData...\n")
sce <- tryCatch({
  zellkonverter::readH5AD(input_path)
}, error = function(e) {
  cat(sprintf("Failed to read AnnData: %s\n", e$message), file = stderr())
  quit(status = 0)
})

# Standardize to logcounts if needed
if (!"logcounts" %in% names(SummarizedExperiment::assays(sce))) {
  cat("Computing logcounts...\n")
  SummarizedExperiment::assay(sce, "logcounts") <- log1p(
    SummarizedExperiment::assay(sce, "counts")
  )
}

# Run SingleR
cat("Running SingleR...\n")
pred <- tryCatch({
  SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main,
    de.method = "classic"
  )
}, error = function(e) {
  cat(sprintf("SingleR prediction failed: %s\n", e$message), file = stderr())
  quit(status = 0)
})

# Add labels to colData
sce$cell_type <- pred$labels
sce$annotation_method <- "SingleR"

cat("Writing output...\n")
zellkonverter::writeH5AD(sce, output_path)

# Write summary
summary <- list(
  method = "singler",
  reference = ref_name,
  cells_labeled = ncol(sce),
  n_unique = length(unique(pred$labels))
)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE), stdout())

cat("SingleR annotation complete.\n")
