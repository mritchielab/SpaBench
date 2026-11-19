#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(sctransform)
  library(patchwork)
  library(ggplot2)
})

option_list <- list(
  make_option("--input",      type="character"),
  make_option("--output",     type="character"),
  make_option("--assay",      type="character", default="Spatial"),
  make_option("--hvg_csv",    type="character"),
  make_option("--plot_png",   type="character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# load, run SCTransform
seu <- readRDS(opt$input)
# ### NEW: VISIUM HD SPECIFIC TISSUE FILTERING
# ==============================================================================
if (!is.null(opt$platform) && opt$platform == "VisiumHD") {
  message("Detected VisiumHD: Filtering for in-tissue bins...")
  
  # Check if slice1 exists and contains the tissue column
  if ("slice1" %in% names(seu@images)) {
    coords <- seu@images$slice1@coordinates
    if ("tissue" %in% colnames(coords)) {
      in_tissue_barcodes <- rownames(coords[coords$tissue == 1, ])
      seu <- subset(seu, cells = in_tissue_barcodes)
      message("Kept ", length(in_tissue_barcodes), " in-tissue bins.")
    } else {
      warning("Tissue column not found in coordinates. Skipping tissue filter.")
    }
  } else {
    warning("Image 'slice1' not found. Skipping tissue filter.")
  }
}

# ### NEW: GLOBAL ZERO-COUNT FILTER
# This prevents log(0) errors (Inf/-Inf) in LogNormalize and SCTransform
seu$nCount_Spatial <- colSums(GetAssayData(seu, assay = "Spatial", slot = "counts"))
if (any(seu$nCount_Spatial == 0)) {
  message("Removing bins with zero total counts to prevent normalization errors...")
  seu <- subset(seu, subset = nCount_Spatial > 0)
}
# ==============================================================================**
seu <- SCTransform(seu, assay = opt$assay, verbose = FALSE)

# export variable features and scaled data
write.csv(VariableFeatures(seu), opt$hvg_csv, row.names = FALSE)
#write.csv(as.data.frame(GetAssayData(seu, slot = "scale.data")), opt$scaled_csv, row.names = TRUE)
saveRDS(seu, opt$output)

# plot top variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot_png, p, width = 14, height = 12)

