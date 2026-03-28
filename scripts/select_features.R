# scripts/select_features.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})
option_list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--method", type="character"),
  make_option("--n_features", type="integer"),
  make_option("--hvg_csv", type="character"),
  make_option("--plot", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
if (!is.null(seu@commands$NormalizeData.Spatial) && seu@commands$NormalizeData.Spatial@params$normalization.method == 'RC'){
  mat <- GetAssayData(seu, assay=assay, slot="data")
  mat <- log1p(mat)
  if ("layers" %in% slotNames(seu[[assay]])) {
    seu[[assay]]@layers$data <- mat
  } else {
    seu[[assay]]$data <- mat
  }
}
seu <- ScaleData(seu, assay=assay, verbose=TRUE)
seu <- tryCatch({
  message("Running FindVariableFeatures on native Assay5...")
  FindVariableFeatures(seu, selection.method = opt$method, nfeatures = opt$n_features)
  
}, error = function(e) {
  # If the above fails, e.g., "argument 1 is not a vector"
  message("Assay5 selection failed. Casting to legacy Assay for compatibility...")
  
  # 3. Cast the specific assay to the v3/v4 'Assay' class then back to 'Assay5' to maintain v5 workflow
  seu[[assay]] <- as(object = seu[[assay]], Class = "Assay")
  seu[[assay]] <- as(object = seu[[assay]], Class = "Assay5")
  
  # 4. Run the selection on the legacy object
  seu <- FindVariableFeatures(seu, selection.method = opt$method, nfeatures = opt$n_features)
  
  return(seu)
})
write.csv(VariableFeatures(seu), opt$hvg_csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot, p)
