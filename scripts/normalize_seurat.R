#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(Matrix)
})

option_list <- list(
  make_option("--input",       type="character", help="Input Seurat RDS"),
  make_option("--output",      type="character", help="Output normalized RDS"),
  make_option("--method",      type="character", help="Normalization method"),
  make_option("--scale_factor",type="numeric",   help="Scale factor"),
  make_option("--plot_png",    type="character", help="PNG file for violin plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$plot_png) || opt$plot_png == "") {
  stop("Please provide --plot_png for normalize_seurat.R")
}

# Load and normalize
raw <- readRDS(opt$input)
#seu <- NormalizeData(raw, 
#                     assay="Spatial",
#                     normalization.method = opt$method,
#                     scale.factor        = opt$scale_factor)
# ==============================================================================
# ### NEW: VISIUM HD SPECIFIC TISSUE FILTERING
# ==============================================================================
if (!is.null(opt$platform) && opt$platform == "VisiumHD") {
  message("Detected VisiumHD: Filtering for in-tissue bins...")
  
  # Check if slice1 exists and contains the tissue column
  if ("slice1" %in% names(raw@images)) {
    coords <- raw@images$slice1@coordinates
    if ("tissue" %in% colnames(coords)) {
      in_tissue_barcodes <- rownames(coords[coords$tissue == 1, ])
      raw <- subset(raw, cells = in_tissue_barcodes)
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
raw$nCount_Spatial <- colSums(GetAssayData(raw, assay = "Spatial", slot = "counts"))
if (any(raw$nCount_Spatial == 0)) {
  message("Removing bins with zero total counts to prevent normalization errors...")
  raw <- subset(raw, subset = nCount_Spatial > 0)
}
# ==============================================================================**
if (opt$method %in% c("LogNormalize", "CLR")) {
  # LogNormalize / CLR 
  seu <- NormalizeData(
    object = raw,
    assay = "Spatial",
    normalization.method = opt$method,
    scale.factor = opt$scale_factor
  )
} else if (opt$method == "RC") {
  # RC 先转 counts 为 dgCMatrix
  cnts <- GetAssayData(raw, assay = "Spatial", slot = "counts")
  if (!inherits(cnts, "dgCMatrix")) {
    cnts <- as(cnts, "dgCMatrix")
    raw <- SetAssayData(
      object = raw,
      assay = "Spatial",
      slot = "counts",
      new.data = cnts
    )
  }
  # 再跑 RC normalization（RC 不使用 scale.factor）
  seu <- NormalizeData(
    object = raw,
    assay = "Spatial",
    normalization.method = "RC"
  )
} else {
  stop("Unsupported normalization method: ", opt$method)
}
# Export counts before/after

# Save normalized object
saveRDS(seu, file = opt$output)

# Violin plots
p1 <- VlnPlot(raw, features = "nCount_Spatial", assay = "Spatial")
p2 <- VlnPlot(seu, features = "nCount_Spatial", assay = "Spatial")
combined <- p1 + p2

ggsave(filename = opt$plot_png,
       plot     = combined,
       width    = 12,
       height   = 10)
