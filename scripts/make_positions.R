#!/usr/bin/env Rscript
# scripts/make_positions.R
#
# Write data/<sample>/01_positions.csv from 01_raw_seurat.rds.
# Format (no header, one row per spot/bead/bin, same order as the Seurat object):
#   barcode, in_tissue, array_row, array_col, pixel_row, pixel_col
# i.e. the 10x tissue_positions_list.csv layout. Consumers:
#   - cluster_SpaGCN.py: cols 3/4 = array (x_array, y_array), cols 5/6 = pixel (x_pixel, y_pixel)
#   - compute_cluster_metrics.R: col 1 = barcode, last two columns = x / y
# Visium / VisiumHD (VisiumV1 image): tissue, row, col, imagerow, imagecol (full-res pixels).
# Slide-seq / Stereo-seq (SlideSeq image): in_tissue = 1 and the bead / bin
#   coordinates (x, y) are used for both the array and the pixel columns.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",  type = "character", help = "01_raw_seurat.rds"),
  make_option("--output", type = "character", help = "01_positions.csv")
)))
if (is.null(opt$input) || is.null(opt$output)) stop("--input and --output are required")

seu <- readRDS(opt$input)
if (length(seu@images) == 0) stop("No spatial image in ", opt$input)
img_name <- if ("slice1" %in% names(seu@images)) "slice1" else names(seu@images)[1]
img <- seu@images[[img_name]]

if (inherits(img, "VisiumV2")) {
  # FOV-based Visium image (Load10X_Spatial): centroids are x = imagerow, y = imagecol
  tc <- GetTissueCoordinates(img, scale = NULL)
  rownames(tc) <- tc$cell
  md <- seu@meta.data
  ar <- if ("array_row" %in% colnames(md)) md[rownames(tc), "array_row"] else tc$x
  ac <- if ("array_col" %in% colnames(md)) md[rownames(tc), "array_col"] else tc$y
  pos <- data.frame(barcode = rownames(tc), in_tissue = 1L,
                    array_row = ar, array_col = ac, pixel_row = tc$x, pixel_col = tc$y)
} else {
  cd <- as.data.frame(img@coordinates)
  cn <- colnames(cd)
  if (all(c("row", "col", "imagerow", "imagecol") %in% cn) && !("x" %in% cn)) {
    # VisiumV1 (Visium / VisiumHD)
    tissue <- if ("tissue" %in% cn) cd$tissue else 1L
    pos <- data.frame(barcode = rownames(cd), in_tissue = tissue,
                      array_row = cd$row, array_col = cd$col,
                      pixel_row = cd$imagerow, pixel_col = cd$imagecol)
  } else {
    # SlideSeq image (Slide-seq beads / Stereo-seq bins)
    x <- if ("x" %in% cn) cd$x else if ("imagecol" %in% cn) cd$imagecol else cd[[1]]
    y <- if ("y" %in% cn) cd$y else if ("imagerow" %in% cn) cd$imagerow else cd[[2]]
    pos <- data.frame(barcode = rownames(cd), in_tissue = 1L,
                      array_row = x, array_col = y, pixel_row = x, pixel_col = y)
  }
}

# keep the Seurat object's cell order; drop coordinates of cells not in the object
pos <- pos[match(intersect(colnames(seu), pos$barcode), pos$barcode), , drop = FALSE]
if (nrow(pos) == 0) stop("No barcodes shared between the image coordinates and the object")
if (anyNA(pos[, 3:6])) stop("Missing coordinates for ", sum(!complete.cases(pos[, 3:6])), " spots")
if (nrow(pos) < ncol(seu))
  warning(ncol(seu) - nrow(pos), " cells of the object have no coordinates and are not written")

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
write.table(pos, opt$output, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Wrote ", nrow(pos), " positions to ", opt$output)
