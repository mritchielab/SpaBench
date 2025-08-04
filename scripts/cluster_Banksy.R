#!/usr/bin/env Rscript
# scripts/cluster_Banksy.R
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(SeuratWrappers)
  library(Banksy)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type="character", help="Input Seurat RDS"),
  make_option("--output",     type="character", help="Output clustered RDS"),
  make_option("--resolution", type="numeric",   help="Clustering resolution"),
  make_option("--meta_csv",   type="character", help="CSV for cluster labels"),
  make_option("--plot_png",   type="character", help="PNG for cluster plot"),
  make_option("--lambda",     type="numeric",   default=0.8, help="Banksy λ"),
  make_option("--k_geom",     type="integer",   default=50,  help="Banksy k_geom")
)))

# 1. load
seu <- readRDS(opt$input)

# 2. figure out feat
# input path .../{feat}/04_reduced.rds
feat <- basename(dirname(opt$input))
message("Feature selection = ", feat)

# 3. set assay
assay <- DefaultAssay(seu)

# 4. ensure var.features for non‑vst SVG methods
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

# 5. run BANKSY
seu <- RunBanksy(
  object   = seu,
  lambda   = opt$lambda,
  verbose  = TRUE,
  assay    = assay,
  slot     = "data",
  features = "variable",
  k_geom   = opt$k_geom
)

# 6. downstream PCA + clustering
DefaultAssay(seu) <- "BANKSY"
seu <- RunPCA(seu, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seu), npcs = 30)
seu <- FindNeighbors(seu, reduction = "pca.banksy", dims = 1:30)
seu <- FindClusters(seu, cluster.name = "banksy_cluster", resolution = opt$resolution)
Idents(seu) <- "banksy_cluster"

# 7. save outputs
write.csv(seu@meta.data["banksy_cluster",drop=FALSE],
          file      = opt$meta_csv,
          row.names = TRUE)
saveRDS(seu, file = opt$output)

# 8. plot
p1 <- SpatialDimPlot(seu, 
                     group.by = "banksy_cluster", 
                     label = T, 
                     repel = T, 
                     label.size = 4, 
                     pt.size.factor = 2)
p2 <- DimPlot(seu, group.by = "banksy_cluster", 
              label = T, 
              repel = T, 
              label.size = 4, 
              pt.size = 2)
ggsave(filename=opt$plot_png,
       plot    = p1 + p2,
       width   = 14,
       height  = 12)
