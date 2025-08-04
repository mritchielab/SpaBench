#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SpatialPCA)
  library(ggplot2)
})

option_list <- list(
  make_option("--input",       type="character", help="Input Seurat RDS"),
  make_option("--output",      type="character", help="Output clustered Seurat RDS"),
  make_option("--meta_csv",    type="character", help="Output CSV of cluster labels"),
  make_option("--plot_png",    type="character", help="PNG for cluster plot"),
  make_option("--n_clusters",  type="integer",   default=7,    help="Number of clusters"),
  make_option("--gene_number", type="integer",   default=2000, help="Number of genes"),
  make_option("--refine",      type="character",        help="If set, run refine step (fast=FALSE)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load Seurat object
seu <- readRDS(opt$input)
feat <- basename(dirname(opt$input))
message("Feature selection = ", feat)
assay <- DefaultAssay(seu)

# 2. Feature‐selection logic (same as BayesSpace example)
if (!(feat %in% c("vst","nnSVG","SPARK","SPARK-X","disp","mvp"))) {
  if (feat == "moransi") {
    message("Overriding VariableFeatures via moransi")
    VariableFeatures(seu, assay=assay) <- SpatiallyVariableFeatures(
      seu, assay=assay, selection.method="moransi"
    )
  } else if (feat == "markvariogram") {
    message("Overriding VariableFeatures via markvariogram")
    VariableFeatures(seu, assay=assay) <- SpatiallyVariableFeatures(
      seu, assay=assay, selection.method="markvariogram"
    )
  } else {
    stop("Unsupported feature‐selection method: ", feat)
  }
} else if (feat %in% c("vst","disp","mvp")) {
  message("Using existing VariableFeatures (vst/disp/mvp)")
} else if (feat %in% c("nnSVG","SPARK","SPARK-X")) {
  message("Using features from ", feat, " (already set)")
} else {
  stop("Unknown feature‐selection method: ", feat)
}

# 3. Prepare spatial coords
coords_df <- as.data.frame(seu@images$slice1@coordinates)[,c("imagecol","imagerow")]
colnames(coords_df) <- c("x_coord","y_coord")
coords <- as.matrix(coords_df)
assay <- DefaultAssay(seu)

# 4. Build SpatialPCA object
SPCA_obj <- CreateSpatialPCAObject(
  counts     = as.matrix(GetAssayData(seu, assay = assay, slot = "counts")),
  location   = coords,
  project    = "SpatialPCA",
  gene.type  = "custom",
  gene.number = opt$gene_number,
  customGenelist = VariableFeatures(seu)
  #min.locations  = 20,
  #min.features   = 20
)

SPCA_obj <- SpatialPCA_buildKernel(
  SPCA_obj,
  kerneltype       = "gaussian",
  bandwidthtype    = "SJ",
  bandwidth.set.by.user = NULL
)

# 5. Estimate loadings & compute Spatial PCs
fast_flag <- if (opt$refine) FALSE else TRUE
SPCA_obj <- SpatialPCA_EstimateLoading(
  SPCA_obj,
  fast         = fast_flag,
  SpatialPCnum = opt$n_clusters * 3  # or use config["pca"]["npcs"] if preferred
)
SPCA_obj <- SpatialPCA_SpatialPCs(
  SPCA_obj,
  fast = fast_flag
)

# 6. Clustering
k_nn <- round(sqrt(ncol(SPCA_obj@SpatialPCs)))
clusterlabel <- walktrap_clustering(
  clusternum = opt$n_clusters,
  latent_dat = SPCA_obj@SpatialPCs,
  knearest   = k_nn
)

# optionally refine
clusterlabel <- refine_cluster_10x(
    clusterlabels = clusterlabel,
    location      = SPCA_obj@location,
    shape         = "hexagon"
  )

# 7. Attach results to Seurat
seu@meta.data$SpatialPCA <- factor(clusterlabel)

# 8. Save outputs
saveRDS(seu, opt$output)
write.csv(
  data.frame(spot=rownames(seu@meta.data), SpatialPCA=clusterlabel),
  file = opt$meta_csv, row.names = FALSE
)

# 9. Plot
#cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91", "chocolate4", "darkorange", "darkcyan", "darkolivegreen2", "lightpink3", "seagreen1", "mediumpurple4")
#p <- plot_cluster(
#  location     = coords,
#  clusterlabel = clusterlabel,
#  pointsize    = 3.5,
#  color_in     = cbp,
#  title_in     = "SpatialPCA"
#)
p <- SpatialDimPlot(
       seu,
       group.by        = "SpatialPCA",
       label           = TRUE,
       repel           = TRUE,
       pt.size.factor  = 2
   )
ggsave(filename = opt$plot_png, plot = p, width=16, height=14)

