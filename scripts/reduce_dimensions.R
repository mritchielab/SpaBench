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
  make_option(c("--umap_dims"), type="integer", default=NULL,
              help="PCs used for UMAP (same as the Seurat clustering); default = all computed PCs"),
  make_option(c("--umap_csv"), type="character", help="Output CSV file for UMAP embeddings"),
  make_option(c("--platform"), type="character", default="Visium", help="Platform: Visium, VisiumHD, Slideseq, Stereoseq"),
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
assay <- DefaultAssay(seu)
if (feat %in% c("moransi", "markvariogram")) {
  message("Using SpatiallyVariableFeatures(method = '", feat, "') for PCA")
  pca_features <- SpatiallyVariableFeatures(seu, method = feat, assay = assay)
} else {
  message("Using VariableFeatures(seu) from '", feat, "' for PCA")
  pca_features <- VariableFeatures(seu, assay = assay)
}
if (is.null(pca_features) || length(pca_features) == 0) {
  stop("No selected features found in the object; please check that feature selection was run successfully.")
}

# 3b. PCA input (scale.data) for the selected features
if (inherits(seu[[assay]], "SCTAssay")) {
  # SCTransform: PCA on the SCT Pearson residuals. The feature-selection step may
  # have overwritten scale.data with ScaleData() (z-scored log1p corrected counts),
  # and SCTransform stores residuals only for its own 3,000 HVGs, so recompute the
  # residuals of all selected features from the SCT model (replace.value = TRUE).
  umi_assay <- slot(slot(seu[[assay]], "SCTModel.list")[[1]], "umi.assay")
  message("SCT assay: computing Pearson residuals for ", length(pca_features),
          " selected features (umi.assay = ", umi_assay, ")")
  seu <- GetResidual(seu, features = pca_features, assay = assay, umi.assay = umi_assay,
                     replace.value = TRUE, verbose = FALSE)
  missing <- setdiff(pca_features, rownames(GetAssayData(seu, assay = assay, layer = "scale.data")))
  if (length(missing) > 0) {
    warning(length(missing), " selected features have no SCT residuals (not in the SCT model) and are excluded from PCA")
    pca_features <- setdiff(pca_features, missing)
  }
} else {
  # Other normalisations: z-score the selected features of the normalised data
  seu <- ScaleData(seu, features = pca_features, verbose = FALSE)
}
seu <- RunPCA(seu, features = pca_features, npcs = opt$npcs)

# 4. Perform UMAP dimensionality reduction on the computed PCs
umap_dims <- min(if (is.null(opt$umap_dims)) opt$npcs else opt$umap_dims, ncol(Embeddings(seu, "pca")))
message("PCA: ", ncol(Embeddings(seu, "pca")), " PCs computed; UMAP uses PCs 1:", umap_dims)
platform <- tolower(opt$platform)
message("Platform: ", opt$platform)

if (platform == "visiumhd") {
  message("VisiumHD detected: using uwot::umap with random init to avoid spectral segfault")
  pca_emb <- Embeddings(seu, "pca")[, seq_len(umap_dims), drop = FALSE]
  
  set.seed(42)
  umap_emb <- uwot::umap(
    pca_emb,
    n_neighbors  = 30L,
    n_components = 2,
    metric       = "cosine",
    init         = "random",
    n_threads    = 8L
  )
  
  rownames(umap_emb) <- rownames(pca_emb)
  colnames(umap_emb) <- c("umap_1", "umap_2")
  seu[["umap"]] <- CreateDimReducObject(
    embeddings = umap_emb,
    key        = "umap_",
    assay      = DefaultAssay(seu)
  )
} else {
  message("Non-VisiumHD platform: using standard RunUMAP")
  seu <- RunUMAP(seu, dims = 1:umap_dims)
}

# 5. Write UMAP embeddings to CSV
write.csv(Embeddings(seu, "umap"), file = opt$umap_csv, row.names = TRUE)

# 6. Slim the object before saving: drop scale.data to keep RDS size manageable
#    (scale.data is only needed for PCA which is already done)
seu <- DietSeurat(seu,
                  counts     = TRUE,
                  data       = TRUE,
                  scale.data = FALSE,
                  dimreducs  = c("pca", "umap"),
                  graphs     = NULL)
saveRDS(seu, file = opt$output)

# 7. Generate and save a UMAP plot
p <- DimPlot(seu, reduction = "umap") 
ggsave(filename = opt$plot_png, plot = p, width = 14, height = 12)
