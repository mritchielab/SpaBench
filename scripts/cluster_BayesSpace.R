#!/usr/bin/env Rscript
# scripts/cluster_BayesSpace.R

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(BayesSpace)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

option_list <- list(
  make_option("--input",      type="character", help="Input Seurat RDS"),
  make_option("--output",     type="character", help="Output clustered Seurat RDS"),
  make_option("--sce_output", type="character", help="Output SingleCellExperiment RDS"),
  make_option("--meta_csv",   type="character", help="CSV for cluster labels"),
  make_option("--plot_png",   type="character", help="PNG for cluster plot"),
  make_option("--nrep",        type="integer",   default=1000,     help="Number of MCMC reps"),
  make_option("--gamma",       type="double",    default=2,        help="Spatial smoothing parameter"),
  make_option("--burn_in",     type="integer",   default=100,      help="MCMC burn-in steps"),
  make_option("--q",           type="integer",   default=10,       help="Number of spatial clusters"),
  make_option("--d",           type="integer",   default=15,       help="Top PCs fed into MCMC (BayesSpace recommends 7-15)"),
  make_option("--model",       type="character", default="t",      help="Error model: 't' (robust, handles outlier spots) or 'normal'"),
  make_option("--init_method", type="character", default="mclust", help="MCMC initialisation: 'mclust' (GMM, better) or 'kmeans' (faster)"),
  make_option("--platform",    type="character", default="Visium", help="Visium, VisiumHD, Stereoseq or Slideseq"),
  make_option("--precision",   type="character", default="equal",  help="Covariance structure: 'equal' (shared) or 'variable' (per-cluster)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load Seurat object
seu <- readRDS(opt$input)

# input path .../{feat}/04_reduced.rds
feat <- basename(dirname(opt$input))
message("Feature selection = ", feat)
assay <- DefaultAssay(seu)

# features used for clustering
feat <- basename(dirname(opt$input))
message("Feature selection = ", feat)
assay <- DefaultAssay(seu)

if (!(feat %in% c("vst", "nnSVG", "SPARK", "SPARK-X", "disp", "mvp"))) {
  # moransi / markvariogram
  if (feat == "moransi") {
    message("Overriding var.features with SpatiallyVariableFeatures(method = 'moransi')")
    VariableFeatures(seu, assay = assay) <- SpatiallyVariableFeatures(
      seu, assay = assay, selection.method = "moransi"
    )
  } else if (feat == "markvariogram") {
    message("Overriding var.features with SpatiallyVariableFeatures(method = 'markvariogram')")
    VariableFeatures(seu, assay = assay) <- SpatiallyVariableFeatures(
      seu, assay = assay, selection.method = "markvariogram"
    )
  } else {
    stop("Unsupported feature‐selection method: ", feat)
  }
  
} else if (feat %in% c("vst", "disp", "mvp")) {
  message("Using existing VariableFeatures (vst / disp / mvp)")
  
} else if (feat %in% c("nnSVG", "SPARK", "SPARK-X")) {
  message(sprintf(
    "Using features from %s (already set in VariableFeatures(seu))",
    feat
  ))
  
} else {
  stop("Unknown feature‐selection method: ", feat)
}
# seu@assays$Spatial@var.features <- SpatiallyVariableFeatures(seu, assay = assay, selection.method=feat)

# 2. Convert to SCE
diet.seurat <- DietSeurat(seu, assays=assay, dimreducs="pca")
sce <- as.SingleCellExperiment(diet.seurat)

# 3. Add spatial coordinates into colData(sce)
# BayesSpace finds neighbours from INTEGER array indices (array_row/array_col).
# The previous code passed pixel coordinates (imagerow/imagecol), so BayesSpace
# found 0 neighbours and ran without its spatial prior.
cd  <- as.data.frame(seu@images[[1]]@coordinates)
cd  <- cd[colnames(sce), , drop = FALSE]
plat <- tolower(opt$platform)
to_grid <- function(v) {                 # integer grid index from evenly spaced coordinates
  u <- sort(unique(v)); step <- min(diff(u)[diff(u) > 0])
  as.integer(round((v - min(v)) / step))
}
if (plat == "visium" && all(c("row", "col") %in% colnames(cd))) {
  coords <- data.frame(array_row = cd$row, array_col = cd$col, row.names = rownames(cd))
  bs_platform <- "Visium"                 # hexagonal grid
} else {
  # VisiumHD bins and Stereo-seq bins lie on a square grid -> BayesSpace "ST" (4 neighbours)
  rc <- if (all(c("row", "col") %in% colnames(cd))) cd[, c("row", "col")] else cd[, c("imagerow", "imagecol")]
  coords <- data.frame(array_row = to_grid(rc[[1]]), array_col = to_grid(rc[[2]]), row.names = rownames(cd))
  bs_platform <- "ST"
  if (plat == "slideseq") warning("Slide-seq beads are not on a grid; BayesSpace neighbours will be sparse.")
}
colData(sce) <- cbind(colData(sce), coords)


# 4. BayesSpace preprocess & cluster
sce <- spatialPreprocess(
  sce,
  platform     = bs_platform,
  skip.PCA     = TRUE,
  log.normalize= FALSE
)
sce <- spatialCluster(
  sce,
  q           = opt$q,
  d           = opt$d,
  gamma       = opt$gamma,
  model       = opt$model,
  init.method = opt$init_method,
  precision   = opt$precision,
  nrep        = opt$nrep,
  burn.in     = opt$burn_in
)

#Add BayesSpace results to Seurat
seu@meta.data = cbind(seu@meta.data, BayesSpace = sce$spatial.cluster)

# 5. Save SCE object
saveRDS(sce, file = opt$sce_output)

# 7. Write cluster labels to CSV
write.csv(seu@meta.data,
          file      = opt$meta_csv,
          row.names = TRUE)

# 8. Save updated Seurat object
saveRDS(seu, file = opt$output)

# 9. Generate and save plots (wrapped to avoid plot failures blocking metrics output)
tryCatch({
  p1 <- SpatialDimPlot(
    seu,
    group.by        = "BayesSpace",
    label           = TRUE,
    repel           = TRUE,
    pt.size.factor  = 2
  )
  ggsave(
    filename = opt$plot_png,
    plot     = p1,
    width    = 14,
    height   = 12
  )
}, error = function(e) {
  message("Plot step failed (non-fatal): ", conditionMessage(e))
  png(opt$plot_png, width = 800, height = 600)
  plot.new(); text(0.5, 0.5, paste("Plot failed:", conditionMessage(e)), cex = 0.8)
  dev.off()
})
