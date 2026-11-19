# scripts/select_features.R
#!/usr/bin/env Rscript
# Highly variable feature selection with Seurat's mean.var.plot (mvp) or
# dispersion (disp) statistics.
#
#  * Both statistics are computed explicitly from the assay's data layer
#    (Seurat:::DISP / CalcDispersion), so the named method is applied for every
#    assay type, including SCTAssay (where FindVariableFeatures() would silently
#    dispatch to the SCT method and ignore selection.method).
#  * mvp : genes ranked by scaled dispersion (mvp.dispersion.scaled), top n_features.
#  * disp: genes ranked by dispersion (mvp.dispersion), top n_features.
#  * The statistics assume log-scale input. For RC (non-log) data a log1p copy
#    is used for ranking only; the stored data layer is never modified.
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

method <- opt$method
if (method %in% c("mean.var.plot", "mvp")) {
  method <- "mean.var.plot"
} else if (method %in% c("dispersion", "disp")) {
  method <- "dispersion"
} else {
  stop("Unsupported --method '", opt$method, "' (expected mean.var.plot or dispersion)")
}

seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
message("Assay: ", assay, " (class ", class(seu[[assay]])[1], "); method: ", method)

# Local copy of the data layer used only for ranking
mat <- LayerData(seu, assay = assay, layer = "data")
if (!inherits(mat, "dgCMatrix")) mat <- as(as.matrix(mat), "CsparseMatrix")
is_rc <- !is.null(seu@commands$NormalizeData.Spatial) &&
  identical(seu@commands$NormalizeData.Spatial@params$normalization.method, "RC")
if (is_rc) {
  message("RC normalisation detected: ranking on a log1p copy of the data layer (stored data unchanged)")
  mat <- log1p(mat)
}

hvf <- Seurat:::CalcDispersion(object = mat, verbose = FALSE)
score_col <- if (method == "mean.var.plot") "mvp.dispersion.scaled" else "mvp.dispersion"
score <- hvf[[score_col]]
names(score) <- rownames(hvf)
expressed <- Matrix::rowSums(mat != 0) > 0
score <- score[expressed[names(score)] & is.finite(score)]
selected <- head(names(sort(score, decreasing = TRUE)), opt$n_features)
rm(mat)
message("Selected ", length(selected), " features (requested ", opt$n_features,
        ", rankable ", length(score), ")")

VariableFeatures(seu[[assay]]) <- selected
stopifnot(identical(VariableFeatures(seu), selected))

write.csv(selected, opt$hvg_csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot mean vs dispersion, highlighting selected features
hvf$selected <- rownames(hvf) %in% selected
p <- ggplot(hvf, aes(x = mvp.mean, y = .data[[score_col]], colour = selected)) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "red")) +
  labs(x = "Average expression (mvp.mean)", y = score_col,
       title = sprintf("%s: %d selected features (%s assay)", method, length(selected), assay)) +
  theme_bw()
ggsave(opt$plot, p)
