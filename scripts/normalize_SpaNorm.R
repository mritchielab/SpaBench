#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(SpaNorm)
  library(SpatialExperiment)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--input", type="character", help="Input Seurat RDS"),
  make_option("--output", type="character", help="Output normalized Seurat RDS"),
  make_option("--spe_output", type="character", help="Output SpatialExperiment RDS"),
  make_option("--plot_png", type="character", help="PNG file for violin plots"),
  make_option("--platform", type="character"),
  make_option("--sample_id", type="character", default="sample1", help="sample id to put into SPE (optional)")
)
opt <- parse_args(OptionParser(option_list = option_list))
platform <- opt$platform

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$spe_output) ||  is.null(opt$plot_png)) {
  stop("Please provide --input, --output, --spe_output, --counts_csv and --plot_png")
}
# Load Seurat object
seu <- readRDS(opt$input)

# transform seu to spe object 
seurat_to_spe <- function(seu, sample_id, img_id) {
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## 1. Build a DataFrame for rowData, pulling in feature names
  feature_names <- rownames(sce)
  rd <- DataFrame(
    gene_id   = feature_names,
    gene_name = feature_names
  )
  rownames(rd) <- feature_names
  
  ## 2. Spatial coords & image as before
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")]
  )
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image)
  )
  imgData <- DataFrame(
    sample_id   = sample_id,
    image_id    = img_id,
    data        = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres
  )
  
  ## 3. Create the SpatialExperiment, injecting your rowData
  spe <- SpatialExperiment(
    assays       = assays(sce),
    rowData      = rd,            # <<< use your new DataFrame here
    colData      = colData(sce),
    metadata     = metadata(sce),
    reducedDims  = reducedDims(sce),
    altExps      = altExps(sce),
    sample_id    = sample_id,
    spatialCoords= spatialCoords,
    imgData      = imgData
  )
  
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  spe
}

seurat_to_spe_stereoslideseq <- function(seu, sample_id, img_id) {
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## 1. Build a DataFrame for rowData, pulling in feature names
  feature_names <- rownames(sce)
  rd <- DataFrame(
    gene_id   = feature_names,
    gene_name = feature_names
  )
  rownames(rd) <- feature_names
  
  ## 2. Spatial coords & image as before
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")]
  )
  
  ## 3. Create the SpatialExperiment, injecting your rowData
  spe <- SpatialExperiment(
    assays       = assays(sce),
    rowData      = rd,            # <<< use your new DataFrame here
    colData      = colData(sce),
    metadata     = metadata(sce),
    reducedDims  = reducedDims(sce),
    altExps      = altExps(sce),
    sample_id    = sample_id,
    spatialCoords= spatialCoords
  )
  
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  spe
}

if (platform %in% c("Visium", "VisiumHD")) {
  spe <- seurat_to_spe(
    seu,           
    sample_id = "data1",
    img_id    = "slice1"
  )
} else if (platform %in% c("Slideseq", "Stereoseq")) {
  spe <- seurat_to_spe_stereoslideseq(
    seu,           
    sample_id = "data1",
    img_id    = "slice1"
  )
  SummarizedExperiment::assay(spe, "counts") <- round(SummarizedExperiment::assay(spe, "counts"))
}

# SpaNorm workflow
keep = filterGenes(spe, 0.2)
table(keep)
spe = spe[keep, ]
# Filter out spots that have zero total counts across the remaining genes
keep_spots <- colSums(counts(spe)) > 0
message("Removing ", sum(!keep_spots), " spots with zero counts.")
spe <- spe[, keep_spots]
logcounts(spe) = log2(counts(spe) + 1)
spe <- SpaNorm(spe)
fit.spanorm <- metadata(spe)$SpaNorm
# Save SpatialExperiment object
saveRDS(spe, file = opt$spe_output)

norm_mat <- as.matrix(assay(spe, "logcounts"))  # genes x cells
orig_counts <- GetAssayData(seu, assay = "Spatial", slot = "counts")
genes <- intersect(rownames(norm_mat), rownames(seu))
cells <- intersect(colnames(norm_mat), colnames(seu))
if (length(genes) == 0 || length(cells) == 0) {
  stop("No overlapping genes or cells between SPE result and Seurat object")
}
norm_sub <- norm_mat[genes, cells, drop = FALSE]
counts_sub <- orig_counts[genes, cells, drop = FALSE]
# backup original Spatial assay
seu@misc$Spatial_backup <- seu[["Spatial"]]
# replace Spatial assay（counts + data）
seu[["Spatial"]] <- CreateAssay5Object(
  counts = counts_sub,
  data   = as.matrix(norm_sub)
)
# write normalized counts CSV (data slot)
#write.csv(as.data.frame(GetAssayData(seu, assay = "Spatial", slot = "data")),
#          file = opt$counts_csv,
#          row.names = TRUE)

# save normalized Seurat object
saveRDS(seu, file = opt$output)

# Violin plots: nCount_Spatial before and after
p1 <- tryCatch(VlnPlot(readRDS(opt$input), features = "nCount_Spatial", assay = "Spatial") , error = function(e) NULL)
p2 <- tryCatch(VlnPlot(readRDS(opt$output), features = "nCount_Spatial", assay = "Spatial") , error = function(e) NULL)

# if Violin plots couldn't be made due to missing metadata, try simple nCount_Spatial violin if exists
if (is.null(p1) || is.null(p2)) {
  # try to compute nCount_Spatial if not present
  obj_raw <- readRDS(opt$input)
  if (!"nCount_Spatial" %in% colnames(obj_raw@meta.data)) {
    obj_raw@meta.data$nCount_Spatial <- Matrix::colSums(GetAssayData(obj_raw, assay = "Spatial", slot = "counts"))
  }
  obj_norm <- readRDS(opt$output)
  if (!"nCount_Spatial" %in% colnames(obj_norm@meta.data)) {
    obj_norm@meta.data$nCount_Spatial <- Matrix::colSums(GetAssayData(obj_norm, assay = "Spatial", slot = "counts"))
  }
  p1 <- VlnPlot(obj_raw, features = "nCount_Spatial", assay = "Spatial")
  p2 <- VlnPlot(obj_norm, features = "nCount_Spatial", assay = "Spatial")
}

combined <- p1 + p2
ggsave(filename = opt$plot_png,
       plot = combined,
       width = 16,
       height = 14)
