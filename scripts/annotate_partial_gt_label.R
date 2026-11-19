#!/usr/bin/env Rscript
#
# annotate_partial_gt_label.R
#
# Annotation + metrics for samples whose ground truth labels only cover a
# subset of spots (e.g. Mouse_brain_VisiumHD: ~15k of ~600k barcodes labelled).
# Both the annotation and the metric calculation are restricted to the
# truth-labelled subset only — unlabelled spots are not annotated and not
# plotted, so the output figure shows just the labelled region.
#
# Pipeline (per truth column):
#   1. Build SingleR reference from labelled cells (per-class average expr).
#   2. Run SingleR on the labelled subset only.
#   3. Per-cluster majority vote of SingleR labels, using labelled cells only.
#   4. Map cluster -> label, assign annotation_<tc> on labelled cells only
#      (unlabelled cells get NA so they drop out of plots).
#   5. Compute NMI / ARI / HOM / COM (annotation vs truth) on labelled cells.
#
# Usage:
#   Rscript scripts/annotate_partial_gt_label.R <sample> <norm> <feat> <alg> [save_rds]
# Example:
#   Rscript scripts/annotate_partial_gt_label.R \
#       Mouse_brain_VisiumHD LogNormalize vst SpatialLeiden
#
# Outputs (under results/<sample>/<norm>/<feat>/<alg>/):
#   09_<alg>_annotated_partial.rds
#   09_<alg>_annotation_partial_<tc>.png        # one per truth column
#   09_<alg>_metrics_partial.csv                # NMI/ARI/HOM/COM per truth col

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(aricode)
  library(mclust)
})

# ---- 1. Args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: annotate_partial_gt_label.R <sample> <norm> <feat> <alg> [save_rds]")
}
sample   <- args[1]
norm     <- args[2]
feat     <- args[3]
alg      <- args[4]
save_rds <- if (length(args) >= 5) tolower(args[5]) %in% c("true", "yes", "1") else TRUE

base_data    <- file.path("data",    sample, norm, feat, alg)
base_results <- file.path("results", sample, norm, feat, alg)

seurat_in    <- file.path(base_data,    sprintf("06_%s_clustered.rds", alg))
clusters_in  <- file.path(base_results, sprintf("07_%s_clusters.csv",  alg))
truth_in     <- file.path("data", "Ground_truth", sample, "labels.tsv")
seurat_out   <- file.path(base_results, sprintf("09_%s_annotated_partial.rds", alg))
plot_out     <- file.path(base_results, sprintf("09_%s_annotation_partial.png", alg))
metrics_out  <- file.path(base_results, sprintf("09_%s_metrics_partial.csv",   alg))

for (p in c(seurat_in, clusters_in, truth_in)) {
  if (!file.exists(p)) stop("Missing input: ", p)
}
dir.create(base_results, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[partial-GT] sample=%s norm=%s feat=%s alg=%s",
                sample, norm, feat, alg))

# ---- 2. Read inputs --------------------------------------------------------
seu <- readRDS(seurat_in)
expr_mat <- as.matrix(GetAssayData(seu, layer = "data"))

clust_df <- read.csv(clusters_in, stringsAsFactors = FALSE)
colnames(clust_df)[1] <- "barcode"
clust_df$barcode <- as.character(clust_df$barcode)

truth_df <- read.delim(truth_in, stringsAsFactors = FALSE)
if ("ID" %in% names(truth_df)) {
  truth_df$barcode <- as.character(truth_df$ID)
} else if ("X" %in% names(truth_df)) {
  truth_df$barcode <- as.character(truth_df$X)
} else {
  stop("Could not find barcode column (ID or X) in truth file.")
}

# ---- 3. Coordinate detection (mirrors annotate_clusters.R) -----------------
image_coords <- NULL
rn <- rownames(seu@meta.data)

if (grepl("^[0-9]+_[0-9]+", rn[1])) {
  message("Platform detected: Stereo-seq. Parsing coordinates from rownames.")
  image_coords <- data.frame(
    col = as.integer(sub("_.*$", "", rn)),
    row = as.integer(sub("^[^_]*_", "", rn)),
    row.names = rn
  )
} else {
  message("Platform detected: Visium/SlideSeq. Checking slots and metadata.")
  coords <- NULL
  try({ coords <- GetTissueCoordinates(seu) }, silent = TRUE)

  if (!is.null(coords) && nrow(coords) > 0) {
    if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
      image_coords <- coords[, c("imagecol", "imagerow")]
    } else {
      image_coords <- coords[, 1:2]
    }
  } else {
    possible_x <- c("x_pixel", "imagecol", "x", "spatial_x", "array_col", "x_array")
    possible_y <- c("y_pixel", "imagerow", "y", "spatial_y", "array_row", "y_array")
    x_col <- intersect(possible_x, colnames(seu@meta.data))[1]
    y_col <- intersect(possible_y, colnames(seu@meta.data))[1]
    if (!is.na(x_col) && !is.na(y_col)) {
      message("Using coordinates from metadata columns: ", x_col, " and ", y_col)
      image_coords <- seu@meta.data[, c(x_col, y_col)]
    }
  }
}

if (is.null(image_coords) || any(is.na(image_coords))) {
  image_coords <- seu@reductions$Xspatial_@cell.embeddings
  colnames(image_coords) <- c("x", "y")
}

# ---- 4. Identify cluster column on the Seurat object ----------------------
if ("refined_pred" %in% colnames(seu@meta.data)) {
  seu@meta.data$cluster <- seu@meta.data$refined_pred
} else if ("spatialleiden" %in% colnames(seu@meta.data)) {
  seu@meta.data$cluster <- seu@meta.data$spatialleiden
} else if ("stlearn" %in% colnames(seu@meta.data)) {
  seu@meta.data$cluster <- seu@meta.data$stlearn
} else {
  banksy_col <- grep("banksy", colnames(seu@meta.data), ignore.case = TRUE, value = TRUE)
  if (length(banksy_col) > 0) {
    seu@meta.data$cluster <- seu@meta.data[[banksy_col[1]]]
  } else {
    seu@meta.data$cluster <- clust_df[[ncol(clust_df)]][
      match(rownames(seu@meta.data), clust_df$barcode)
    ]
  }
}

# ---- 5. Barcode harmonization ---------------------------------------------
if (length(intersect(colnames(seu), truth_df$barcode)) == 0) {
  message("Exact barcode match failed. Attempting to harmonize suffixes...")
  seu_clean   <- sub("[-_][0-9]$|-[0-9]$|_[0-9]$", "", colnames(seu))
  truth_clean <- sub("[-_][0-9]$|-[0-9]$|_[0-9]$", "", truth_df$barcode)
  if (length(intersect(seu_clean, truth_clean)) > 0) {
    seu <- RenameCells(seu, new.names = seu_clean)
    expr_mat <- GetAssayData(seu, layer = "data")
    truth_df$barcode <- truth_clean
    clust_df$barcode <- sub("[-_][0-9]$|-[0-9]$|_[0-9]$", "", clust_df$barcode)
    rownames(image_coords) <- seu_clean
  } else {
    stop("Barcodes do not match even after suffix cleaning.")
  }
}

# ---- 6. Image slot --------------------------------------------------------
colnames(image_coords) <- c("col", "row")
seu[["image"]] <- new(Class = "SlideSeq",
                     assay = "Spatial",
                     coordinates = as.data.frame(image_coords))

# ---- 7. Per-truth-column annotation (labelled subset only) ---------------
truth_cols  <- setdiff(names(truth_df), c("X", "ID", "barcode"))
cluster_col <- tail(names(clust_df), 1)
metrics_rows <- list()

for (tc in truth_cols) {
  message("Processing truth column: ", tc)

  label_map <- truth_df[[tc]]
  names(label_map) <- truth_df$barcode
  groups <- label_map[colnames(seu)]

  valid_cells <- !is.na(groups) & nzchar(groups)
  n_labeled   <- sum(valid_cells)
  message(sprintf("  labelled cells: %d / %d", n_labeled, length(groups)))
  if (n_labeled == 0) next

  # Reference: per-class average expression over labelled cells only.
  sub_expr   <- expr_mat[, valid_cells, drop = FALSE]
  sub_groups <- groups[valid_cells]
  group_indices <- split(seq_len(ncol(sub_expr)), sub_groups)
  ref_avg <- do.call(cbind, lapply(group_indices, function(idx) {
    rowMeans(sub_expr[, idx, drop = FALSE])
  }))

  se_ref <- SummarizedExperiment(
    assays  = SimpleList(logcounts = as.matrix(ref_avg)),
    colData = DataFrame(label = colnames(ref_avg), row.names = colnames(ref_avg))
  )

  # SingleR runs on the labelled subset only — we never assign labels outside it.
  sr <- SingleR(test = sub_expr, ref = se_ref, labels = se_ref$label,
                assay.type.test = "logcounts", assay.type.ref = "logcounts")

  pred_df <- data.frame(barcode  = colnames(sub_expr),
                        sr_label = sr$labels,
                        stringsAsFactors = FALSE)

  # Per-cluster majority vote across labelled cells only.
  merged <- clust_df %>% inner_join(pred_df, by = "barcode")
  if (nrow(merged) == 0) next

  annot <- merged %>%
    group_by(cluster = !!sym(cluster_col)) %>%
    count(sr_label) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup()

  # Apply cluster->label map ONLY to labelled cells; rest stays NA so plots
  # and metrics ignore them.
  meta_col <- paste0("annotation_", tc)
  cluster_vec <- seu@meta.data$cluster
  full_pred <- annot$sr_label[match(cluster_vec, annot$cluster)]
  full_pred[!valid_cells] <- NA
  seu@meta.data[[meta_col]] <- factor(full_pred)

  # ---- Metrics on labelled subset ----------------------------------------
  true_v <- as.character(groups[valid_cells])
  pred_v <- as.character(full_pred[valid_cells])
  ok     <- !is.na(pred_v)
  true_v <- true_v[ok]; pred_v <- pred_v[ok]

  if (length(true_v) >= 2 &&
      length(unique(true_v)) >= 2 &&
      length(unique(pred_v)) >= 2) {
    homogeneity  <- function(t, p) { x <- table(p, t); sum(apply(x, 1, max)) / length(t) }
    completeness <- function(t, p) { x <- table(p, t); sum(apply(x, 2, max)) / length(t) }
    metrics_rows[[tc]] <- data.frame(
      truth     = tc,
      n_labeled = n_labeled,
      n_scored  = length(true_v),
      NMI       = aricode::NMI(true_v, pred_v),
      ARI       = mclust::adjustedRandIndex(true_v, pred_v),
      HOM       = homogeneity(true_v, pred_v),
      COM       = completeness(true_v, pred_v),
      stringsAsFactors = FALSE
    )
  } else {
    message("  skipping metrics (single class after filtering)")
  }
}

# ---- 8. Save annotated RDS ------------------------------------------------
if (save_rds) {
  saveRDS(seu, file = seurat_out)
  message("Saved annotated Seurat object to: ", seurat_out)
} else {
  writeLines("Annotation RDS saving skipped.", seurat_out)
  message("Wrote placeholder to: ", seurat_out)
}

# ---- 9. Metrics CSV -------------------------------------------------------
if (length(metrics_rows) > 0) {
  metrics_df <- do.call(rbind, metrics_rows)
  write.csv(metrics_df, metrics_out, row.names = FALSE)
  message("Wrote metrics CSV: ", metrics_out)
} else {
  message("No metrics computed — nothing written.")
}

# ---- 10. Plots: only labelled cells (subset Seurat) -----------------------
for (tc in truth_cols) {
  meta_col <- paste0("annotation_", tc)
  if (!(meta_col %in% colnames(seu@meta.data))) next
  labelled_cells <- colnames(seu)[!is.na(seu@meta.data[[meta_col]])]
  if (length(labelled_cells) < 2) next

  seu_sub <- subset(seu, cells = labelled_cells)
  seu_sub@meta.data[[meta_col]] <- droplevels(factor(seu_sub@meta.data[[meta_col]]))

  p <- SpatialDimPlot(seu_sub, group.by = meta_col, label = TRUE) +
    labs(title = paste("Cell Type Annotation (partial GT, Ref:", tc, ")"))
  plot_path <- sub("\\.png$", paste0("_", tc, ".png"), plot_out)
  ggsave(plot_path, p, width = 16, height = 14)
}

message("Annotation (partial GT) completed.")
