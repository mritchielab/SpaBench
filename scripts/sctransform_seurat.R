#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(sctransform)
  library(patchwork)
  library(ggplot2)
})

option_list <- list(
  make_option("--input",      type="character"),
  make_option("--output",     type="character"),
  make_option("--assay",      type="character", default="Spatial"),
  make_option("--hvg_csv",    type="character"),
  make_option("--scaled_csv", type="character"),
  make_option("--plot_png",   type="character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# load, run SCTransform
seu <- readRDS(opt$input)
seu <- SCTransform(seu, assay = opt$assay, verbose = FALSE)

# export variable features and scaled data
write.csv(VariableFeatures(seu), opt$hvg_csv, row.names = FALSE)
write.csv(as.data.frame(GetAssayData(seu, slot = "scale.data")), opt$scaled_csv, row.names = TRUE)
saveRDS(seu, opt$output)

# plot top variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot_png, p, width = 14, height = 12)

