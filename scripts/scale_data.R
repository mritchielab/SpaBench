# scripts/scale_data.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})
option_list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--scaled_csv", type="character"),
  make_option("--plot_png", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

seu <- readRDS(opt$input)
seu <- ScaleData(seu, features=VariableFeatures(seu))
write.csv(as.data.frame(GetAssayData(seu, slot="scale.data")), opt$scaled_csv, row.names=TRUE)
saveRDS(seu, opt$output)

# plot scaled feature distribution
p <- VlnPlot(seu, features=VariableFeatures(seu)[1:5])
ggsave(opt$plot_png, p)