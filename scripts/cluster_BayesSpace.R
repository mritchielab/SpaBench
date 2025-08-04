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
  make_option("--nrep",       type="integer",   default=1000, help="Number of MCMC reps"),
  make_option("--burn_in",    type="integer",   default=100,  help="MCMC burn-in"),
  make_option("--q",          type="integer",   default=10,   help="Number of clusters")
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
diet.seurat <- DietSeurat(seu, assays=assay, graphs="pca")
sce <- as.SingleCellExperiment(diet.seurat)

# 3. Add spatial coordinates into colData(sce)
# rename 'imagerow' as 'array_row'
coords <- seu@images$slice1@coordinates %>%
  dplyr::rename(array_row=imagerow, array_col=imagecol)
colData(sce) <- cbind(colData(sce), coords)


# 4. BayesSpace preprocess & cluster
sce <- spatialPreprocess(
  sce,
  platform     = "Visium",
  skip.PCA     = TRUE,
  log.normalize= FALSE
)
sce <- spatialCluster(
  sce,
  q       = opt$q,
  #d       = NULL,
  #hist.max= NULL,
  #gamma   = NULL,
  #init.method = "mclust",
  #model   = "t", 
  nrep    = opt$nrep,
  burn.in = opt$burn_in
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

# 9. Generate and save plots
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
