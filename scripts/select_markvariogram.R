#!/usr/bin/env Rscript
# Spatially variable feature selection with the mark variogram
# (Seurat::FindSpatiallyVariableFeatures, selection.method = "markvariogram").
#
#  * Candidates: all genes detected (data layer > 0) in at least 3 spots, NOT a
#    vst pre-selection. Optionally capped with --n_candidates (most frequently
#    detected genes kept) for very large datasets.
#  * Candidates are ranked by the mark-variogram r.metric (Seurat ranking) and the top n_features kept.
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

sv_method <- "markvariogram"

# Define command-line options
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type="character", help="Path to input Seurat RDS file"),
  make_option("--output",     type="character", help="Path to save updated Seurat RDS file"),
  make_option("--csv",        type="character", help="Path to output selected feature CSV"),
  make_option("--method",     type="character", help="Set as markvariogram"),
  make_option("--plot",       type="character", help="Path to save feature plot PNG"),
  make_option("--n_features", type="integer",   default=2000,
              help="Number of spatially variable features to keep"),
  make_option("--n_candidates", type="integer", default=0,
              help="Max number of candidate genes to rank (0 = all genes detected in >= 3 spots); if capped, the most frequently detected genes are kept"),
  make_option("--workers", type="integer", default=1,
              help="Parallel worker processes (forked) for the mark-variogram computation (1 = sequential)"),
  make_option("--chunk_size", type="integer", default=100,
              help="Genes per parallel mark-variogram task (bounds per-worker memory)")
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

# Optional parallel pre-computation. Seurat caches per-gene mark variograms in
# Misc(assay, "markvariogram") and FindSpatiallyVariableFeatures() skips genes already
# present there, so we compute them in small forked chunks (bounded memory, unlike
# Seurat's future backend which gives each worker 1/N of all genes at once) and let
# Seurat compute the remaining chunk, merge, rank and flag as usual.
if (!is.na(opt$workers) && opt$workers > 1) {
  img <- SeuratObject:::DefaultImage(seu)  # same image FindSpatiallyVariableFeatures uses
  tc <- GetTissueCoordinates(seu[[img]])
  mv_data <- as.matrix(LayerData(seu[[assay]], layer = stat_slot,
                                 cells = rownames(tc), features = candidates))
  mv_data <- mv_data[Seurat:::RowVar(mv_data) > 0, , drop = FALSE]
  genes <- rownames(mv_data)
  chunk_size <- max(2L, opt$chunk_size)
  chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
  if (length(chunks) > 1) {
    last <- chunks[[length(chunks)]]
    if (length(last) < 2) {  # Seurat's own call needs >= 2 genes left
      chunks[[length(chunks) - 1]] <- head(chunks[[length(chunks) - 1]], -1)
    }
    pre <- chunks[-length(chunks)]
    message("Pre-computing mark variograms for ", length(unlist(pre)), " genes in ",
            length(pre), " chunks on ", opt$workers, " workers")
    mv <- parallel::mclapply(pre, function(g) {
      Seurat:::RunMarkVario(spatial.location = tc, data = mv_data[g, , drop = FALSE])
    }, mc.cores = opt$workers, mc.preschedule = FALSE)
    failed <- which(vapply(seq_along(mv), function(i) {
      x <- mv[[i]]
      is.null(x) || inherits(x, "try-error") || length(x) != length(pre[[i]])
    }, logical(1)))
    if (length(failed)) {
      # e.g. a forked worker killed for memory: recompute those chunks sequentially
      message("Retrying ", length(failed), " failed chunk(s) sequentially")
      for (i in failed) {
        mv[[i]] <- Seurat:::RunMarkVario(spatial.location = tc,
                                         data = mv_data[pre[[i]], , drop = FALSE])
      }
    }
    mv <- unlist(mv, recursive = FALSE)
    names(mv) <- unlist(pre)
    Misc(seu[[assay]], slot = "markvariogram") <- mv
    rm(mv)
  }
  rm(mv_data)
}

# Select spatially variable features using the mark variogram
# (pass `slot`: FindSpatiallyVariableFeatures.Seurat does not forward `layer`)
seu <- FindSpatiallyVariableFeatures(
  object = seu,
  assay = assay,
  slot = stat_slot,
  selection.method = sv_method,
  features = candidates,
  nfeatures = opt$n_features
)
# drop the cached per-gene variogram objects (large, not used downstream)
suppressWarnings(Misc(seu[[assay]], slot = "markvariogram") <- NULL)

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
