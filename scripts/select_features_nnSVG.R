# scripts/select_features_nnSVG.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(nnSVG)
  library(SpatialExperiment)
  library(dplyr)
  library(scran)
})
option_list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--n_features", type="integer"),
  make_option("--hvg_csv", type="character"),
  make_option("--platform", type="character"),
  make_option("--plot_png", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
platform <- opt$platform

seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
seu <- ScaleData(seu, assay=assay, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=opt$n_features)

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
}

# use SpatialExperiment object for nnSVG workflow
filtered_spe <- filter_genes(spe)
if (length(filtered_spe@rowRanges) >= 2000 ) {
  spe <- filtered_spe
}
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
assayNames(spe)
# set seed for reproducibility
set.seed(123)
spe <- nnSVG(spe, n_threads = 1)
table(rowData(spe)$padj <= 0.05)
rowData(spe)[rowData(spe)$padj <= 0.05, ]

# add SVG back to seu object
# Keep significant genes in nnSVG rank order (not genome order), at most n_features, no NA padding
rd <- rowData(spe)
rd <- rd[!is.na(rd$padj) & rd$padj <= 0.05, , drop = FALSE]
VariableFeatures(seu) <- head(rownames(rd)[order(rd$rank)], opt$n_features)
write.csv(VariableFeatures(seu), opt$hvg_csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot_png, p)



