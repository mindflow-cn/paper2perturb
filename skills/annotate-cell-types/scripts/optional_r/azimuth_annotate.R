#!/usr/bin/env Rscript
#
# Azimuth reference-based cell-type mapping (optional adapter).
#
# Usage:
#   Rscript azimuth_annotate.R input.h5ad output.h5ad --reference skin
#
# Reference options: pbmc, skin, lung, kidney, pancreas, liver, brain, etc.
#
# Only invoked when --mode deep --enable-r-adapters are both set.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript azimuth_annotate.R input.h5ad output.h5ad [--reference reference_name]")
}

input_path   <- args[1]
output_path  <- args[2]
ref_name     <- if (length(args) >= 4 && args[3] == "--reference") args[4] else "pbmc"

cat(sprintf("Azimuth annotation: input=%s, output=%s, ref=%s\n",
            input_path, output_path, ref_name))

required_packages <- c("Azimuth", "Seurat", "zellkonverter", "reticulate")
missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing) > 0) {
  msg <- sprintf(
    "Azimuth adapter skipped: missing R packages: %s. Install with:\n  install.packages(c(%s))",
    paste(missing, collapse = ", "),
    paste(sprintf("'%s'", missing), collapse = ", ")
  )
  cat(msg, "\n", file = stderr())
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
})

# Load AnnData and convert to Seurat
cat("Loading AnnData...\n")
obj <- tryCatch({
  zellkonverter::readH5AD(input_path)
}, error = function(e) {
  cat(sprintf("Failed to read AnnData: %s\n", e$message), file = stderr())
  quit(status = 0)
})

# Azimuth reference URL mapping
ref_urls <- list(
  pbmc    = "https://azimuth.hub.consortium.io/reference/pbmc_multimodal.h5seurat",
  skin    = "https://azimuth.hub.consortium.io/reference/skin.h5seurat",
  lung    = "https://azimuth.hub.consortium.io/reference/lung.h5seurat",
  kidney  = "https://azimuth.hub.consortium.io/reference/kidney.h5seurat",
  pancreas = "https://azimuth.hub.consortium.io/reference/pancreas.h5seurat"
)

ref_url <- ref_urls[[ref_name]]
if (is.null(ref_url)) {
  cat(sprintf("Unknown Azimuth reference: %s. Available: %s\n",
              ref_name, paste(names(ref_urls), collapse = ", ")),
      file = stderr())
  quit(status = 0)
}

cat(sprintf("Running Azimuth with reference: %s\n", ref_name))
pred <- tryCatch({
  Azimuth::RunAzimuth(obj, reference = ref_url)
}, error = function(e) {
  cat(sprintf("Azimuth prediction failed: %s\n", e$message), file = stderr())
  quit(status = 0)
})

cat("Writing output...\n")
zellkonverter::writeH5AD(pred, output_path)

summary <- list(
  method = "azimuth",
  reference = ref_name,
  cells_labeled = ncol(pred)
)
if ("jsonlite" %in% installed.packages()[,"Package"]) {
  writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE), stdout())
}

cat("Azimuth annotation complete.\n")
