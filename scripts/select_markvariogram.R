#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

orig_f <- getFromNamespace("FindSpatiallyVariableFeatures.StdAssay", "Seurat")
b <- body(orig_f)
idx <- which(sapply(seq_along(b), function(i) {
  expr_lines <- deparse(b[[i]])
  any(grepl("object\\[names\\(x = svf.info\\)\\]", expr_lines))
}))
if (length(idx) != 1) {
  stop("Check duplicated code!")
}
b[[idx]] <- quote(
  object[[ names(x = svf.info) ]] <- svf.info
)
body(orig_f) <- b
unlockBinding("FindSpatiallyVariableFeatures.StdAssay", asNamespace("Seurat"))
assign("FindSpatiallyVariableFeatures.StdAssay", orig_f, asNamespace("Seurat"))
lockBinding("FindSpatiallyVariableFeatures.StdAssay", asNamespace("Seurat"))

# Define command‑line options
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type = "character", help = "Path to input Seurat RDS file"),
  make_option("--output",     type = "character", help = "Path to save updated Seurat RDS file"),
  make_option("--method", type="character"),
  make_option("--csv",        type = "character", help = "Path to output selected feature CSV"),
  make_option("--plot",       type = "character", help = "Path to save feature plot PNG"),
  make_option("--n_features", type = "integer",   default = NA,
              help = "(Placeholder) Number of features, accepted but not used by this script")
)))

# Load Seurat object
seu <- readRDS(opt$input)

# Ensure we're operating on the 'Spatial' assay, then scale the data
assay <- DefaultAssay(seu)
seu <- ScaleData(seu, assay=assay, verbose=FALSE)

# Identify spatially variable features using mark‑variogram
sv <- FindSpatiallyVariableFeatures(
  object = seu,
  assay = assay,
  slot = "scale.data",
  selection.method = "markvariogram",
  features = VariableFeatures(seu)[1:opt$n_features]
)

# Write out feature list and updated object
write.csv(SpatiallyVariableFeatures(sv, method = "markvariogram"),   file = opt$csv, row.names = FALSE)
saveRDS(sv,    file = opt$output)

# Plot top 10 features
p <- SpatialFeaturePlot(sv, features = head(SpatiallyVariableFeatures(sv, method = "markvariogram"), 10))
ggsave(filename = opt$plot, plot = p, width = 14, height = 12)
