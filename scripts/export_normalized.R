#!/usr/bin/env Rscript
# Export the normalised expression of a Step-1 object (02_<norm>_normalized.rds)
# so that the Python clustering methods (SpaGCN, stLearn, SpatialLeiden) use the
# SAME normalised values as the R methods of that workflow.
#
# Output: HDF5 file with the "data" layer of the object's default assay, stored
# as a compressed-sparse-column matrix (genes x spots):
#   /data /indices /indptr  CSC components (0-based indices)
#   /shape                  c(n_genes, n_spots)
#   /genes /barcodes        row / column names
#   /assay                  assay name the values came from
# Read it with scripts/norm_data.py.
suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(Matrix); library(hdf5r)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",  type = "character", help = "02_<norm>_normalized.rds"),
  make_option("--output", type = "character", help = "Output .h5")
)))

seu   <- readRDS(opt$input)
assay <- DefaultAssay(seu)
m     <- GetAssayData(seu, assay = assay, layer = "data")
m     <- as(m, "CsparseMatrix")
if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
message(sprintf("Exporting assay '%s' data layer: %d genes x %d spots", assay, nrow(m), ncol(m)))

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
if (file.exists(opt$output)) file.remove(opt$output)
f <- H5File$new(opt$output, mode = "w")
f[["data"]]     <- as.numeric(m@x)
f[["indices"]]  <- as.integer(m@i)
f[["indptr"]]   <- as.integer(m@p)
f[["shape"]]    <- as.integer(dim(m))
f[["genes"]]    <- rownames(m)
f[["barcodes"]] <- colnames(m)
f[["assay"]]    <- assay
f$close_all()
message("Wrote ", opt$output)
