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
  make_option("--n_clusters",  type="integer",   default=20,    help="Number of clusters"),
  make_option("--gene_number", type="integer",   default=2000, help="Number of genes"),
  make_option("--spatial_pcs", type="integer",   default=20,   help="Number of spatial PCs (SpatialPCnum; package default 20)"),
  make_option("--refine",      type="character",        help="If TRUE, use fast=FALSE loadings/PCs and apply refine_cluster_10x; if FALSE, use fast=TRUE and skip refine_cluster_10x")
)
opt <- parse_args(OptionParser(option_list=option_list))
do_refine <- isTRUE(as.logical(opt$refine))

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

# 4. Build SpatialPCA object from the WORKFLOW's normalised values.
# CreateSpatialPCAObject() would re-normalise the raw counts with its own
# SCTransform (so the Step-1 normalisation would have no effect) and drop spots
# with < 20 detected genes. Instead we build the object directly and set
# normalized_expr to the row-scaled 'data' layer of the default assay
# (= Step-1 normalisation) for the selected genes, which is exactly the slot
# SpatialPCA_buildKernel / EstimateLoading read (SpatialPCA 1.3.0 source:
# CreateSpatialPCAObject stores centred SCT scale.data there and buildKernel
# row-scales object@normalized_expr into object@params$expr).
genes <- VariableFeatures(seu, assay = assay)
norm_data <- GetAssayData(seu, assay = assay, layer = "data")
genes <- intersect(as.character(genes), rownames(norm_data))
if (length(genes) == 0) stop("None of the selected features are in the '", assay, "' data layer")
expr <- as.matrix(norm_data[genes, , drop = FALSE])
rm(norm_data)
# zero-variance genes cannot be scaled (NaN) -> drop them
gene_sd <- apply(expr, 1, sd)
if (any(!is.finite(gene_sd) | gene_sd == 0)) {
  message(sprintf("Dropping %d constant genes", sum(!is.finite(gene_sd) | gene_sd == 0)))
  expr <- expr[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
}
expr <- t(scale(t(expr)))
spots <- intersect(colnames(expr), rownames(coords))
expr <- expr[, spots, drop = FALSE]
coords <- coords[spots, , drop = FALSE]
message(sprintf("SpatialPCA input: '%s' assay data layer (normalisation of this workflow), %d genes x %d spots",
                assay, nrow(expr), ncol(expr)))

SPCA_obj <- new(Class = "SpatialPCA",
                counts   = as(GetAssayData(seu, assay = assay, layer = "counts")[rownames(expr), spots, drop = FALSE], "dgCMatrix"),
                location = coords,
                project  = "SpatialPCA")
SPCA_obj@normalized_expr <- expr
SPCA_obj@params <- list()
rm(expr)

SPCA_obj <- SpatialPCA_buildKernel(
  SPCA_obj,
  kerneltype       = "gaussian",
  bandwidthtype    = "SJ",
  bandwidth.set.by.user = NULL
)

# 5. Estimate loadings & compute Spatial PCs
SPCA_obj <- SpatialPCA_EstimateLoading(
  SPCA_obj,
  fast         = !do_refine,
  SpatialPCnum = opt$spatial_pcs
)
message("SpatialPCA: ", opt$spatial_pcs, " spatial PCs")
SPCA_obj <- SpatialPCA_SpatialPCs(
  SPCA_obj,
  fast = !do_refine
)

# 6. Clustering
k_nn <- round(sqrt(ncol(SPCA_obj@SpatialPCs)))
clusterlabel <- walktrap_clustering(
  clusternum = opt$n_clusters,
  latent_dat = SPCA_obj@SpatialPCs,
  knearest   = k_nn
)

# optionally refine
if (do_refine) {
  clusterlabel <- refine_cluster_10x(
    clusterlabels = clusterlabel,
    location      = SPCA_obj@location,
    shape         = "hexagon"
  )
}

# 7. Attach results to Seurat
# SpatialPCs columns follow colnames(normalized_expr) (= spots)
cl <- setNames(as.character(clusterlabel), spots)
seu@meta.data$SpatialPCA <- factor(cl[rownames(seu@meta.data)])

# 8. Save outputs
saveRDS(seu, opt$output)
write.csv(
  data.frame(spot=spots, SpatialPCA=clusterlabel),
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

