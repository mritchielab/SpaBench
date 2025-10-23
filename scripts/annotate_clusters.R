## scripts/annotate_clusters.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--seurat",     type="character", help="Input Seurat RDS"),
  make_option("--clusters",   type="character", help="Clustering CSV"),
  make_option("--truth",      type="character", help="Ground-truth TSV"),
  make_option("--seurat_out", type="character", help="Output annotated Seurat RDS"),
  make_option("--plot_seurat",type="character", help="Output spatial annotation PNG")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Read Seurat object and data matrix
seu      <- readRDS(opt$seurat)
expr_mat <- as.matrix(GetAssayData(seu, slot="data"))

# 2. Load cluster labels
clust_df <- read.csv(opt$clusters, stringsAsFactors=FALSE)
colnames(clust_df)[1] <- "barcode"
if (any(c("pred", "refined_pred") %in% colnames(seu@meta.data))) {
  seu@meta.data$cluster <- seu@meta.data$refined_pred
  image_coords <- seu@meta.data[, c("x_pixel","y_pixel"), drop = FALSE]
  colnames(image_coords) <- c("col", "row")
  seu[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = image_coords)
} else if (any(c("stlearn") %in% colnames(seu@meta.data))) {
  seu@meta.data$cluster <- seu@meta.data$stlearn
  # if both imagecol and imagerow exist, take them directly
  if (all(c("imagecol", "imagerow") %in% colnames(seu@meta.data))) {
    image_coords <- seu@meta.data[, c("imagecol", "imagerow"), drop = FALSE]
  } else {
    # fallback: parse rownames like "109_198" -> imagecol=109, imagerow=198
    rn <- rownames(seu@meta.data)
    if (is.null(rn)) stop("rownames(seu@meta.data) is NULL, cannot infer image coordinates")
    
    # extract part before first '_' and part after first '_' (safe if multiple underscores)
    imagecol <- as.integer(sub("_.*$", "", rn))
    imagerow  <- as.integer(sub("^[^_]*_", "", rn))
    
    # if parsing failed for many rows, warn the user
    if (all(is.na(imagecol)) || all(is.na(imagerow))) {
      warning("Failed to parse imagecol/imagerow from rownames; parsed values contain mostly NAs")
    }
    
    # assemble as data.frame and also add to meta.data for later use
    image_coords <- data.frame(imagecol = imagecol, imagerow = imagerow, row.names = rn, stringsAsFactors = FALSE)
    seu@meta.data$imagecol <- image_coords$imagecol
    seu@meta.data$imagerow  <- image_coords$imagerow
  }
  seu[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = image_coords)
} else if ("SpatialPCA" %in% colnames(seu@meta.data)) {
  seu@meta.data$cluster <- seu@meta.data$SpatialPCA
  coords0 <- seu@images[[1]]@coordinates
  image_coords <- coords0[rownames(seu@meta.data), c("row","col"), drop = FALSE]
  colnames(image_coords) <- c("col","row")
  seu[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = image_coords)
} else if ("spatialleiden" %in% colnames(seu@meta.data)) {
  seu@meta.data$cluster <- seu@meta.data$spatialleiden
  # if both imagecol and imagerow exist, take them directly
  if (all(c("array_col","array_row") %in% colnames(seu@meta.data))) {
    image_coords <- seu@meta.data[, c("array_col","array_row"), drop = FALSE]
  } else {
    # fallback: parse rownames like "109_198" -> imagecol=109, imagerow=198
    rn <- rownames(seu@meta.data)
    if (is.null(rn)) stop("rownames(seu@meta.data) is NULL, cannot infer image coordinates")
    
    # extract part before first '_' and part after first '_' (safe if multiple underscores)
    imagecol <- as.integer(sub("_.*$", "", rn))
    imagerow  <- as.integer(sub("^[^_]*_", "", rn))
    
    # if parsing failed for many rows, warn the user
    if (all(is.na(imagecol)) || all(is.na(imagerow))) {
      warning("Failed to parse imagecol/imagerow from rownames; parsed values contain mostly NAs")
    }
    
    # assemble as data.frame and also add to meta.data for later use
    image_coords <- data.frame(imagecol = imagecol, imagerow = imagerow, row.names = rn, stringsAsFactors = FALSE)
    seu@meta.data$imagecol <- image_coords$imagecol
    seu@meta.data$imagerow  <- image_coords$imagerow
  }
  seu[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = image_coords)
}

# 3. Load ground truth and build reference
truth_df <- read.delim(opt$truth, stringsAsFactors=FALSE)
if ("ID" %in% names(truth_df)) {
  truth_df$barcode <- truth_df$ID
} else if ("X" %in% names(truth_df)) {
  truth_df$barcode <- truth_df$X
} else {
  stop("Cannot find ID or X column in truth file")
}

truth_cols <- setdiff(names(truth_df), c("X","ID","barcode"))
cluster_col <- tail(names(clust_df),1)
cluster_maps <- list()

for(tc in truth_cols) {
  # build pseudo-bulk reference
  groups  <- truth_df[[tc]]
  ref_avg <- do.call(cbind, lapply(split(seq_len(ncol(expr_mat)), groups),
                      function(idx) rowMeans(expr_mat[, idx, drop=FALSE])))
  colnames(ref_avg) <- names(split(seq_len(ncol(expr_mat)), groups))
  rownames(ref_avg) <- rownames(expr_mat)
  se_ref <- SummarizedExperiment(
    assays  = SimpleList(logcounts=ref_avg),
    colData = DataFrame(label=colnames(ref_avg), row.names=colnames(ref_avg))
  )

  # run SingleR
  sr <- SingleR(test=expr_mat, ref=se_ref, labels=se_ref$label,
                assay.type.test="logcounts", assay.type.ref="logcounts")
  pred_df <- data.frame(barcode=colnames(expr_mat), sr_label=sr$labels,
                        stringsAsFactors=FALSE)

  # merge cluster, truth and prediction
  merged <- clust_df %>%
    inner_join(truth_df, by="barcode") %>%
    inner_join(pred_df, by="barcode")

  # map cluster -> majority sr_label
  annot <- merged %>%
    group_by_at(cluster_col) %>%
    count(sr_label) %>%
    slice_max(n, with_ties=FALSE) %>%
    select(!!sym(cluster_col), annot=sr_label)
  cluster_maps[[tc]] <- annot

  # add annotation_<tc> to metadata
  meta_col <- paste0("annotation_", tc)
  seu@meta.data[[meta_col]] <- factor(
    annot$annot[ match(seu@meta.data[[cluster_col]], annot[[cluster_col]]) ]
  )
}

# save annotated Seurat
dir.create(dirname(opt$seurat_out), showWarnings=FALSE, recursive=TRUE)
saveRDS(seu, file=opt$seurat_out)

# spatial plot per truth
library(ggplot2)
library(patchwork)
plots <- list()
for(tc in truth_cols) {
  meta_col <- paste0('annotation_',tc)
  p <- SpatialDimPlot(seu, group.by=meta_col, label=TRUE) + labs(title=paste('Annotation by',tc))
  ggsave(sub('\\.(rds)$', paste0('_',tc,'.png'),opt$plot_seurat), p, width=16, height=14)
}
