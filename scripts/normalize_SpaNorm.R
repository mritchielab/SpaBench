#!/usr/bin/env Rscript
# SpaNorm normalisation.
#
# Like the other normalisations (normalize_seurat.R, sctransform_seurat.R) this
#   * keeps ALL genes (feature selection happens downstream), and
#   * keeps ALL spots except those with zero total counts (removed with subset()
#     on the whole Seurat object, exactly as the other normalisation scripts do).
# The Spatial assay keeps its original counts; only its 'data' layer is filled
# with the SpaNorm (logPAC) normalised values.
#
# SpaNorm's model is fitted on every gene with at least one count. Genes that
# SpaNorm cannot normalise are not dropped:
#   * all-zero genes are 0 (their value under any log normalisation);
#   * genes with counts whose SpaNorm values are constant or non-finite (SpaNorm
#     learnt nothing for them, e.g. no counts in the spots sampled for fitting)
#     get LogNormalize values (scale factor 1e4) instead, so the gene universe
#     is identical to the other normalisations.
#
# The model is fitted with SpaNorm's own fitting routine (identical defaults to
# SpaNorm::SpaNorm(), adj.method = "auto" = logPAC); the logPAC adjustment is then
# applied in gene chunks into a sparse matrix. SpaNorm::SpaNorm() would build
# several dense genes x spots matrices, which is infeasible for VisiumHD with all
# genes. The chunked adjustment reproduces SpaNorm:::normaliseLogPAC exactly
# (each gene's computation only depends on that gene, except the dispersion cap
# psi.max, which is computed once over all genes as SpaNorm does).
suppressPackageStartupMessages({
  library(optparse)
  library(SpaNorm)
  library(SpatialExperiment)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--input", type="character", help="Input Seurat RDS"),
  make_option("--output", type="character", help="Output normalized Seurat RDS"),
  make_option("--spe_output", type="character", help="Output SpatialExperiment RDS"),
  make_option("--plot_png", type="character", help="PNG file for violin plots"),
  make_option("--platform", type="character"),
  make_option("--sample_id", type="character", default="sample1", help="sample id to put into SPE (optional)"),
  make_option("--ncores",    type="integer",   default=1L,         help="cores (unused: SpaNorm 1.2 fits single-threaded)")
)
opt <- parse_args(OptionParser(option_list = option_list))
platform <- opt$platform

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$spe_output) || is.null(opt$plot_png)) {
  stop("Please provide --input, --output, --spe_output and --plot_png")
}

seu <- readRDS(opt$input)
img_id <- "slice1"

# ---------------------------------------------------------------------------
# 1. Spot filter: identical to normalize_seurat.R / sctransform_seurat.R
# ---------------------------------------------------------------------------
seu$nCount_Spatial <- colSums(GetAssayData(seu, assay = "Spatial", layer = "counts"))
p1 <- VlnPlot(seu, features = "nCount_Spatial", assay = "Spatial")
if (any(seu$nCount_Spatial == 0)) {
  message("Removing ", sum(seu$nCount_Spatial == 0), " spots with zero total counts.")
  seu <- subset(seu, subset = nCount_Spatial > 0)
}

counts <- GetAssayData(seu, assay = "Spatial", layer = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

# SpaNorm's negative binomial model requires integer counts (Slide-seq/Stereo-seq
# matrices may be stored as doubles); round for fitting only.
fit_counts <- counts
fit_counts@x <- round(fit_counts@x)
fit_counts <- drop0(fit_counts)

# A spot whose rounded counts are all zero cannot be handled by SpaNorm
# (log library size = -Inf); remove it from the whole object.
zero_rounded <- colSums(fit_counts) == 0
if (any(zero_rounded)) {
  message("Removing ", sum(zero_rounded), " spots whose rounded counts are all zero (SpaNorm cannot fit them).")
  seu <- subset(seu, cells = colnames(fit_counts)[!zero_rounded])
  counts <- counts[, !zero_rounded, drop = FALSE]
  fit_counts <- fit_counts[, !zero_rounded, drop = FALSE]
}
stopifnot(identical(colnames(counts), Cells(seu)))
message(sprintf("Normalising %d genes x %d spots", nrow(counts), ncol(counts)))

# ---------------------------------------------------------------------------
# 2. SpatialExperiment with the genes that have counts
# ---------------------------------------------------------------------------
fit_genes <- rownames(fit_counts)[Matrix::rowSums(fit_counts) > 0]
message(sprintf("Fitting SpaNorm on %d genes with counts (%d all-zero genes are set to 0).",
                length(fit_genes), nrow(fit_counts) - length(fit_genes)))

coords <- as.matrix(seu@images[[img_id]]@coordinates[colnames(counts), c("imagecol", "imagerow")])
spe_args <- list(
  assays        = list(counts = fit_counts[fit_genes, , drop = FALSE]),
  rowData       = S4Vectors::DataFrame(gene_id = fit_genes, gene_name = fit_genes, row.names = fit_genes),
  colData       = S4Vectors::DataFrame(seu@meta.data[colnames(counts), , drop = FALSE]),
  sample_id     = opt$sample_id,
  spatialCoords = coords
)
if (platform %in% c("Visium", "VisiumHD") && length(seu@images[[img_id]]@image) > 0) {
  spe_args$imgData <- S4Vectors::DataFrame(
    sample_id   = opt$sample_id,
    image_id    = img_id,
    data        = I(list(SpatialExperiment::SpatialImage(as.raster(seu@images[[img_id]]@image)))),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres
  )
}
spe <- do.call(SpatialExperiment, spe_args)
spe$in_tissue <- 1
rm(fit_counts); invisible(gc())

# ---------------------------------------------------------------------------
# 3. Fit SpaNorm (same arguments/defaults as SpaNorm::SpaNorm())
# ---------------------------------------------------------------------------
# For VisiumHD use a small sample.p - 5% of ~500k spots still gives ~25k for robust fitting
sample_p <- if (platform == "VisiumHD") 0.05 else 0.25
message(sprintf("Fitting SpaNorm %s: sample.p=%.2f", packageVersion("SpaNorm"), sample_p))
SpaNorm:::checkSPE(spe)
set.seed(123)   # SpaNorm samples spots for fitting; seed for reproducibility
fit <- SpaNorm:::fitSpaNorm(
  Y = assay(spe, "counts"), coords = spatialCoords(spe), sample.p = sample_p,
  gene.model = "nb", msgfun = message, df.tps = 6L, lambda.a = 1e-4, batch = NULL,
  LS = SingleCellExperiment::sizeFactors(spe), tol = 1e-4, step.factor = 0.5,
  maxit.nb = 50, maxit.psi = 25, maxn.psi = 500
)
metadata(spe)$SpaNorm <- fit

# ---------------------------------------------------------------------------
# 4. logPAC adjustment (SpaNorm:::normaliseLogPAC, scale.factor = 1), gene-chunked
# ---------------------------------------------------------------------------
scale_factor <- 1
psi <- fit$psi
psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))   # global cap, as in normaliseLogPAC
psi <- pmin(psi, psi.max)
isbio <- fit$wtype %in% "biology"
Y <- assay(spe, "counts")
chunk_size <- max(1L, floor(5e7 / ncol(Y)))
chunks <- split(seq_len(nrow(Y)), ceiling(seq_len(nrow(Y)) / chunk_size))
message(sprintf("Applying logPAC adjustment in %d gene chunks", length(chunks)))
norm_fit <- do.call(rbind, lapply(chunks, function(idx) {
  mu   <- SpaNorm:::calculateMu(fit$gmean[idx], fit$alpha[idx, , drop = FALSE], fit$W)
  mu.2 <- SpaNorm:::calculateMu(fit$gmean[idx], fit$alpha[idx, isbio, drop = FALSE], fit$W[, isbio])
  y    <- as.matrix(Y[idx, , drop = FALSE])
  size <- 1 / psi[idx]
  lb <- pnbinom(y - 1, mu = mu, size = size)
  ub <- dnbinom(y, mu = mu, size = size) + lb
  p  <- pmax(pmin((lb + ub) / 2, 0.999), 0.001)
  out <- log2(qnbinom(p, mu = scale_factor * mu.2, size = size) + 1)
  out[!is.finite(out)] <- NA_real_   # flagged below and replaced
  as(out, "CsparseMatrix")
}))
dimnames(norm_fit) <- dimnames(Y)

# Genes SpaNorm could not normalise: non-finite values, or constant output
# (all of the gene's information was lost) -> LogNormalize values instead.
gene_na  <- Matrix::rowSums(is.na(norm_fit)) > 0
rng      <- sparseMatrixStats::rowRanges(norm_fit, na.rm = TRUE)
fallback <- fit_genes[gene_na | rng[, 1] == rng[, 2]]
message(sprintf("%d of %d fitted genes could not be normalised by SpaNorm (non-finite or constant output); using LogNormalize values for them.",
                length(fallback), length(fit_genes)))

# ---------------------------------------------------------------------------
# 5. Assemble the full gene x spot data matrix (sparse) in the assay's order
# ---------------------------------------------------------------------------
ok_genes   <- setdiff(fit_genes, fallback)
zero_genes <- setdiff(rownames(counts), fit_genes)
parts <- list(norm_fit[ok_genes, , drop = FALSE])
if (length(fallback) > 0) {
  ln <- counts[fallback, , drop = FALSE] %*% Diagonal(x = 1e4 / colSums(counts))  # LogNormalize
  ln@x <- log1p(ln@x)
  parts <- c(parts, list(as(ln, "CsparseMatrix")))
}
if (length(zero_genes) > 0) {
  parts <- c(parts, list(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                              dims = c(length(zero_genes), ncol(counts)),
                                              dimnames = list(zero_genes, colnames(counts)))))
}
norm_full <- do.call(rbind, parts)[rownames(counts), , drop = FALSE]
norm_full <- drop0(as(norm_full, "dgCMatrix"))
stopifnot(!anyNA(norm_full@x), all(is.finite(norm_full@x)),
          identical(dimnames(norm_full), dimnames(counts)))
rm(norm_fit, parts); invisible(gc())

# SpatialExperiment output (genes with counts; logcounts = final values used in Seurat)
assay(spe, "logcounts") <- norm_full[fit_genes, , drop = FALSE]
rowData(spe)$SpaNorm_fallback_LogNormalize <- fit_genes %in% fallback
saveRDS(spe, file = opt$spe_output)
rm(spe); invisible(gc())

# Seurat output: original counts (all genes, same spots), SpaNorm values in 'data'
seu <- SetAssayData(seu, assay = "Spatial", layer = "data", new.data = norm_full)
seu@misc$SpaNorm_fallback_LogNormalize <- fallback
stopifnot(identical(Cells(seu[["Spatial"]]), Cells(seu)),
          identical(rownames(seu[["Spatial"]]), rownames(counts)))
saveRDS(seu, file = opt$output)

# Violin plots: nCount_Spatial before / after the spot filter
p2 <- VlnPlot(seu, features = "nCount_Spatial", assay = "Spatial")
ggsave(filename = opt$plot_png, plot = p1 + p2, width = 16, height = 14)
