# scripts/select_features_SPARK-X.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(SPARK)
  library(ggplot2)
  library(Seurat)
})

option_list <- list(
  make_option("--input",      type="character", help="Path to input Seurat RDS file"),
  make_option("--output",     type="character", help="Path to save updated Seurat RDS file"),
  make_option("--n_features", type="integer",   help="Number of top spatial genes to select"),
  make_option("--hvg_csv",    type="character", help="Path to output selected feature CSV"),
  make_option("--plot_png",   type="character", help="Path to save p-value distribution PNG")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Load Seurat object
seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
seu <- ScaleData(seu, assay=assay, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=opt$n_features)

# Extract count matrix and spatial coordinates
sp_count <- as.matrix(GetAssayData(seu, assay = assay, slot = "counts"))
image_name <- names(seu@images)[1]
coords <- as.matrix(seu@images[[image_name]]@coordinates[, c("row", "col")])

# Filter mitochondrial genes if present
mt_idx <- grep("^mt-", rownames(sp_count), ignore.case = TRUE)
if (length(mt_idx) > 0) {
  sp_count <- sp_count[-mt_idx, , drop = FALSE]
}

# Run SPARK-X analysis
data_info <- coords
sparkx_res <- sparkx(sp_count, data_info, numCores = 8, option = "mixture")

# Extract and rank by combined p-value
res_df <- as.data.frame(sparkx_res$res_mtest)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$combinedPval), ]
# Select top n_features
top_genes <- head(res_df$gene, opt$n_features)
# assign SVG selected by SPARK-X to Seurat
VariableFeatures(seu) <- top_genes

write.csv(data.frame(gene = top_genes),
          file   = opt$hvg_csv,
          row.names = FALSE)
saveRDS(seu, file = opt$output)

# 8) p-value distribution
p <- ggplot(res_df, aes(x = -log10(combinedPval))) +
  geom_histogram(bins = 50) +
  labs(title = "SPARK-X Combined P-Value Distribution",
       x = "-log10(Combined P-Value)", y = "Count")
ggsave(opt$plot_png, plot = p, width = 16, height = 14)


