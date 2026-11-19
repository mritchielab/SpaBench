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
  make_option("--input",        type="character",  help="Input Seurat RDS"),
  make_option("--output",       type="character",  help="Output clustered RDS"),
  make_option("--resolution",   type="numeric",    help="Clustering resolution"),
  make_option("--meta_csv",     type="character",  help="CSV for cluster labels"),
  make_option("--plot_png",     type="character",  help="PNG for cluster plot"),
  make_option("--lambda",       type="numeric",    default=0.8,          help="Banksy λ (0=non-spatial, 0.2=cell typing, 0.8=domain)"),
  make_option("--k_geom",       type="integer",    default=50,           help="Banksy k_geom: geometric neighbours"),
  make_option("--use_agf",      action="store_true", default=FALSE,      help="Use Asymmetric Gaussian Function neighbour weighting (pair with --M)"),
  make_option("--M",            type="integer",    default=NULL,         help="Highest azimuthal Fourier harmonic order (only when --use_agf)"),
  make_option("--spatial_mode", type="character",  default="kNN_median", help="Neighbour aggregation kernel: kNN_median, kNN_mean, kNN_r, kNN_rn, rNN_gauss, kNN_unif"),
  make_option("--npcs",         type="integer",    default=20,           help="PCs for RunPCA and FindNeighbors dims"),
  make_option("--k_neighbors",  type="integer",    default=50,           help="sNN graph neighbours for FindNeighbors (k.param)"),
  make_option("--smooth_labels", action="store_true", default=FALSE,     help="Apply spatial k-NN label smoothing after clustering"),
  make_option("--smooth_k",      type="integer",    default=6,           help="k spatial neighbours for label smoothing (default 6)")
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
banksy_args <- list(
  object       = seu,
  lambda       = opt$lambda,
  verbose      = TRUE,
  assay        = assay,
  slot         = "data",
  features     = "variable",
  k_geom       = opt$k_geom,
  use_agf      = opt$use_agf,
  spatial_mode = opt$spatial_mode
)
if (!is.null(opt$M)) banksy_args$M <- opt$M
seu <- do.call(RunBanksy, banksy_args)

# 6. downstream PCA + clustering
DefaultAssay(seu) <- "BANKSY"
seu <- RunPCA(seu, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seu), npcs = opt$npcs)
seu <- FindNeighbors(seu, reduction = "pca.banksy", dims = 1:opt$npcs, k.param = opt$k_neighbors)
seu <- FindClusters(seu, cluster.name = "banksy_cluster", resolution = opt$resolution)
Idents(seu) <- "banksy_cluster"

# 6b. optional spatial label smoothing (majority vote among k spatial neighbours)
# Stored in banksy_cluster_smooth; banksy_cluster left intact for Seurat plotting
if (opt$smooth_labels) {
  suppressPackageStartupMessages(library(FNN))
  img_name <- names(seu@images)[1]
  coords   <- seu@images[[img_name]]@coordinates[, c("imagecol", "imagerow")]
  coords   <- coords[colnames(seu), ]
  nn_idx   <- FNN::get.knnx(as.matrix(coords), as.matrix(coords),
                             k = opt$smooth_k + 1L)$nn.index[, -1L, drop = FALSE]
  raw_lab  <- as.character(seu@meta.data$banksy_cluster)
  smooth_lab <- vapply(seq_len(nrow(nn_idx)), function(i) {
    nb <- raw_lab[nn_idx[i, ]]
    names(sort(table(nb), decreasing = TRUE))[1L]
  }, character(1L))
  seu@meta.data$banksy_cluster_smooth <- smooth_lab
  message("smoothLabels applied: k=", opt$smooth_k,
          "; unique clusters after smoothing = ", length(unique(smooth_lab)))
}

# 7. save outputs — use smoothed labels if available, else raw
out_col <- if (opt$smooth_labels) "banksy_cluster_smooth" else "banksy_cluster"
write.csv(seu@meta.data[out_col],
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
