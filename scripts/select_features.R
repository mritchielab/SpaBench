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
seu <- FindVariableFeatures(seu, selection.method=opt$method, nfeatures=opt$n_features)
write.csv(VariableFeatures(seu), opt$hvg_csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot, p)
