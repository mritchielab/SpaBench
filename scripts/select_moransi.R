#!/usr/bin/env Rscript
# Spatially variable feature selection with Moran's I
# (Seurat::FindSpatiallyVariableFeatures, selection.method = "moransi").
#
#  * Candidates: all genes detected (data layer > 0) in at least 3 spots, NOT a
#    vst pre-selection. Optionally capped with --n_candidates (most frequently
#    detected genes kept) for very large datasets.
#  * Candidates are ranked by Moran's I (Seurat ranking) and the top n_features kept.
#  * Standard assays: statistic computed on scale.data of the candidates (as before).
#    SCTAssay: computed on the SCT data layer (log1p corrected counts).
#  * The selected genes are flagged as spatially variable (SpatiallyVariableFeatures)
#    and also set as VariableFeatures; scale.data is (re)computed for exactly the
#    selected genes, which is what reduce_dimensions.R uses for PCA.
#  * The data layer of the input object is never modified.
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

sv_method <- "moransi"

# Define command-line options
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type="character", help="Path to input Seurat RDS file"),
  make_option("--output",     type="character", help="Path to save updated Seurat RDS file"),
  make_option("--csv",        type="character", help="Path to output selected feature CSV"),
  make_option("--method",     type="character", help="Set as moransi"),
  make_option("--plot",       type="character", help="Path to save feature plot PNG"),
  make_option("--n_features", type="integer",   default=2000,
              help="Number of spatially variable features to keep"),
  make_option("--n_candidates", type="integer", default=0,
              help="Max number of candidate genes to rank (0 = all genes detected in >= 3 spots); if capped, the most frequently detected genes are kept")
)))

# Load Seurat object
seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
is_sct <- inherits(seu[[assay]], "SCTAssay")
message("Assay: ", assay, " (class ", class(seu[[assay]])[1], ")")

# Candidate genes: detected in >= 3 spots
dat <- LayerData(seu, assay = assay, layer = "data")
n_detected <- if (inherits(dat, "dgCMatrix")) {
  Matrix::rowSums(dat > 0)
} else {
  rowSums(as.matrix(dat) > 0)
}
names(n_detected) <- rownames(seu[[assay]])
rm(dat)
candidates <- names(n_detected)[n_detected >= 3]
if (!is.na(opt$n_candidates) && opt$n_candidates > 0 && length(candidates) > opt$n_candidates) {
  candidates <- head(candidates[order(-n_detected[candidates])], opt$n_candidates)
}
message("Ranking ", length(candidates), " candidate genes by ", sv_method)

if (is_sct) {
  stat_slot <- "data"
} else {
  stat_slot <- "scale.data"
  seu <- ScaleData(seu, assay = assay, features = candidates, verbose = FALSE)
}

# Select spatially variable features using Moran's I
# (pass `slot`: FindSpatiallyVariableFeatures.Seurat does not forward `layer`)
seu <- FindSpatiallyVariableFeatures(
  object = seu,
  assay = assay,
  slot = stat_slot,
  selection.method = sv_method,
  features = candidates,
  nfeatures = opt$n_features
)
# NOTE: no x.cuts/y.cuts. With Seurat 5.3 BinData() returns NA positions for these
# coordinates, which made every Moran's I NA (all p-values tied) so the "ranking"
# silently fell back to assay row order. Moran's I is computed per spot instead.
obs <- seu[[assay]][["MoransI_observed"]][, 1]
if (all(is.na(obs))) stop("All Moran's I values are NA; spatial coordinates problem")

selected <- SpatiallyVariableFeatures(seu, method = sv_method)
message("Selected ", length(selected), " features (requested ", opt$n_features, ")")
VariableFeatures(seu[[assay]]) <- selected
# scale.data for exactly the selected genes (used by RunPCA in reduce_dimensions.R)
seu <- ScaleData(seu, assay = assay, features = selected, verbose = FALSE)
stopifnot(identical(SpatiallyVariableFeatures(seu, method = sv_method), selected),
          identical(VariableFeatures(seu), selected))

# Write out feature list and updated object
write.csv(selected, file = opt$csv, row.names = FALSE)
saveRDS(seu, file = opt$output)

# Plot top 10 features
p <- SpatialFeaturePlot(seu, features = head(selected, 10))
ggsave(filename = opt$plot, plot = p, width = 14, height = 12)
