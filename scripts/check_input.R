#!/usr/bin/env Rscript
# check_input.R - minimal version
# Usage: Rscript check_input.R /path/to/xxx.rds

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript check_input.R /path/to/xxx.rds")
inp <- args[[1]]
if (!file.exists(inp)) stop("Input not found: ", inp)
if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat required")

seu <- readRDS(inp)
if (!inherits(seu, "Seurat")) stop("Not a Seurat object")
if (is.null(seu@images) || !("slice1" %in% names(seu@images))) stop("seu@images[['slice1']] missing")

coords <- seu@images[["slice1"]]@coordinates
if (is.null(coords)) stop("coordinates NULL")
coords_df <- as.data.frame(coords, stringsAsFactors = FALSE)
if (!is.null(rownames(coords))) rownames(coords_df) <- rownames(coords)

lc <- tolower(colnames(coords_df))
has_imagepair <- all(c("imagecol","imagerow") %in% lc)
has_colrow   <- all(c("col","row") %in% lc)
has_xy       <- all(c("x","y") %in% lc)

get_col <- function(target) {
  nm <- colnames(coords_df)[which(lc == target)]
  if (length(nm)) nm[1] else NULL
}

if (has_imagepair && has_colrow) {
  cat("input seu 正常\n"); quit(status=0)
}

if (has_imagepair && !has_colrow) {
  imgc <- get_col("imagecol"); imgr <- get_col("imagerow")
  coords_df$col <- coords_df[[imgc]]; coords_df$row <- coords_df[[imgr]]
  msg <- "Copied imagecol->col and imagerow->row"
} else if (has_colrow && !has_imagepair) {
  ccol <- get_col("col"); crow <- get_col("row")
  coords_df$imagecol <- coords_df[[ccol]]; coords_df$imagerow <- coords_df[[crow]]
  msg <- "Copied col->imagecol and row->imagerow"
} else if (!has_imagepair && !has_colrow && has_xy) {
  xn <- get_col("x"); yn <- get_col("y")
  coords_df$imagecol <- coords_df[[xn]]; coords_df$col <- coords_df[[xn]]
  coords_df$imagerow <- coords_df[[yn]]; coords_df$row <- coords_df[[yn]]
  msg <- "Assigned x->(imagecol,col) and y->(imagerow,row)"
} else {
  stop("coordinates lack required columns (imagecol/imagerow or col/row or x/y). Found: ",
       paste(colnames(coords_df), collapse=", "))
}

seu@images[["slice1"]]@coordinates <- coords_df
out <- file.path(dirname(inp), "updated_seu.rds")
saveRDS(seu, out)
cat(msg, "\nSaved updated Seurat to:", out, "\n")

