#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--input"), type="character", help="Path to input Seurat RDS file"),
  make_option(c("--output"), type="character", help="Path to output Seurat RDS file"),
  make_option(c("--npcs"), type="integer", help="Number of principal components to compute"),
  make_option(c("--umap_csv"), type="character", help="Output CSV file for UMAP embeddings"),
  make_option(c("--plot_png"), type="character", help="Output PNG file for UMAP plot")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 1. Load the Seurat object from RDS
seu <- readRDS(opt$input)

# 2. Infer feature selection method from input path
#    Input path format: data/{sample}/{normalization}/{feature_selection}/03_{feature_selection}.rds
feat <- basename(dirname(opt$input))
message("Detected feature selection method: ", feat)

# 3. Determine PCA features based on feature selection method
if (feat == "moransi") {
  message("Using VariableFeatures(seu) from '", feat, "' for PCA")
  pca_features <- SpatiallyVariableFeatures(seu, method = "moransi")
  if (is.null(pca_features) || length(pca_features) == 0) {
    stop("No VariableFeatures found in the object; please check that feature selection was run successfully.")
  }
  # Run PCA using the selected variable features
  seu <- RunPCA(seu, features = pca_features, npcs = opt$npcs)
} else if (feat == "markvariogram") {
  message("Using VariableFeatures(seu) from '", feat, "' for PCA")
  print(DefaultAssay(seu))
  pca_features <- SpatiallyVariableFeatures(seu, method = "markvariogram", assay = DefaultAssay(seu))
  if (is.null(pca_features) || length(pca_features) == 0) {
    stop("No VariableFeatures found in the object; please check that feature selection was run successfully.")
  }
  # Run PCA using the selected variable features
  seu <- RunPCA(seu, features = pca_features, npcs = opt$npcs)
} else {
  message("Feature selection is 'vst'; using default PCA settings")
  # Run PCA with default feature set (Seurat's default behavior)
  assay <- DefaultAssay(seu)
  seu <- ScaleData(seu, assay=assay, verbose=FALSE)
  seu <- RunPCA(seu, npcs = opt$npcs)
}

# 4. Perform UMAP dimensionality reduction on the computed PCs
seu <- RunUMAP(seu, dims = 1:opt$npcs)

# 5. Write UMAP embeddings to CSV
write.csv(Embeddings(seu, "umap"), file = opt$umap_csv, row.names = TRUE)

# 6. Save the updated Seurat object with PCA/UMAP results
saveRDS(seu, file = opt$output)

# 7. Generate and save a UMAP plot
p <- DimPlot(seu, reduction = "umap") 
ggsave(filename = opt$plot_png, plot = p, width = 14, height = 12)
