# scripts/select_features_SPARK.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(SPARK)
})

option_list <- list(
  make_option("--input",      type="character", help="Path to input Seurat RDS file"),
  make_option("--output",     type="character", help="Path to save updated Seurat RDS file"),
  make_option("--n_features", type="integer",   help="Number of top spatial genes to select"),
  make_option("--hvg_csv",    type="character", help="Path to output selected feature CSV"),
  make_option("--plot_png",   type="character", help="Path to save p-value distribution PNG")
)
opt <- parse_args(OptionParser(option_list=option_list))

seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
seu <- ScaleData(seu, assay=assay, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=opt$n_features)
counts   <- GetAssayData(seu, assay = assay, slot = "counts")

# get image name
image_name <- names(seu@images)[1]
coords_df  <- seu@images[[image_name]]@coordinates[, c("row", "col")]
# SPARK metadata
info <- data.frame(x = coords_df[ , "row"],
                   y = coords_df[ , "col"],
                   total_counts = colSums(counts),
                   row.names = colnames(counts))

# 3) create SPARKObject and filter spot
spark <- CreateSPARKObject(counts      = counts,
                           location    = info[, 1:2],
                           percentage  = 0.1,
                           min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
stopifnot(length(spark@lib_size) == ncol(spark@counts))

# 4) fit null model & spatial test
spark <- spark.vc(spark,
                  covariates = NULL,
                  lib_size   = spark@lib_size,
                  num_core   = 5,
                  verbose    = FALSE)
spark <- spark.test(spark,
                    check_positive = TRUE,
                    verbose        = FALSE)

# 5) rank based on adjusted p value as top n_features
res <- spark@res_mtest
res <- res[order(res[ , "adjusted_pvalue"]), , drop = FALSE]
top_genes <- rownames(res)[1:opt$n_features]

# 6) assign SVG selected by SPARK to Seurat
VariableFeatures(seu) <- top_genes

write.csv(data.frame(gene = top_genes),
          file   = opt$hvg_csv,
          row.names = FALSE)
saveRDS(seu, file = opt$output)

# 8) p-value distribution
p <- ggplot(as.data.frame(res), aes(x = adjusted_pvalue)) +
  geom_histogram(binwidth = 0.01) +
  labs(title = "SPARK adjusted p-value distribution",
       x = "adjusted p-value",
       y = "count")
ggsave(opt$plot_png, plot = p, width = 16, height = 14)
