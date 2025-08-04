# scripts/select_vst.R
#!/usr/bin/env Rscript
suppressMessages(library(optparse)); library(Seurat); library(ggplot2)
option_list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--method", type="character"),
  make_option("--n_features", type="integer"),
  make_option("--csv", type="character"),
  make_option("--plot", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
seu <- ScaleData(seu, assay=assay, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=opt$n_features)
write.csv(VariableFeatures(seu), opt$csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot variable features
p <- VariableFeaturePlot(seu)
ggsave(opt$plot, p)
