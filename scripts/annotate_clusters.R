#!/usr/bin/env Rscript

# 1. Load Libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# 2. Parse Arguments
option_list <- list(
  make_option("--seurat",      type="character", help="Input Seurat RDS"),
  make_option("--clusters",    type="character", help="Clustering CSV"),
  make_option("--truth",       type="character", help="Ground-truth TSV"),
  make_option("--seurat_out",  type="character", help="Output annotated Seurat RDS"),
  make_option("--save_rds", type="character", default="True", help="Whether to save annotated RDS"),
  make_option("--plot_seurat", type="character", help="Output spatial annotation PNG")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 3. Read Seurat and Prep Assay Data
seu <- readRDS(opt$seurat)

# Use 'layer' instead of 'slot' for Seurat v5 compatibility
expr_mat <- as.matrix(GetAssayData(seu, layer="data"))

# 4. Load cluster labels
clust_df <- read.csv(opt$clusters, stringsAsFactors=FALSE)
colnames(clust_df)[1] <- "barcode"
clust_df$barcode <- as.character(clust_df$barcode)

# --- MULTI-PLATFORM COORDINATE LOGIC ---
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
      message(paste("Using coordinates from metadata columns:", x_col, "and", y_col))
      image_coords <- seu@meta.data[, c(x_col, y_col)]
    }
  }
}

if (is.null(image_coords) || any(is.na(image_coords))) {
  image_coords <- seu@reductions$Xspatial_@cell.embeddings
  colnames(image_coords) <- c("x", "y")
  #stop("CRITICAL ERROR: No valid spatial coordinates found for this platform.")
}

# --- CLUSTER COLUMN IDENTIFICATION ---
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
    seu@meta.data$cluster <- clust_df[[ncol(clust_df)]][match(rownames(seu@meta.data), clust_df$barcode)]
  }
}

# 5. Load ground truth and Harmonize Barcodes
truth_df <- read.delim(opt$truth, stringsAsFactors=FALSE)
if ("ID" %in% names(truth_df)) {
  truth_df$barcode <- as.character(truth_df$ID)
} else if ("X" %in% names(truth_df)) {
  truth_df$barcode <- as.character(truth_df$X)
} else {
  stop("Error: Could not find barcode column (ID or X) in truth file.")
}

# --- BARCODE HARMONIZATION & COORDINATE SYNC ---
if (length(intersect(colnames(seu), truth_df$barcode)) == 0) {
  message("Exact barcode match failed. Attempting to harmonize suffixes...")
  
  # Strip -1 or _1 common in Visium/Squidpy
  seu_clean <- sub("-[0-9]+$", "", colnames(seu))
  truth_clean <- sub("-[0-9]+$", "", truth_df$barcode)
  
  if (length(intersect(seu_clean, truth_clean)) > 0) {
    # 1. Update Seurat cell names
    seu <- RenameCells(seu, new.names = seu_clean)
    
    # 2. Update expression matrix names
    expr_mat <- GetAssayData(seu, layer="data") 
    
    # 3. Update internal dataframes
    truth_df$barcode <- truth_clean
    clust_df$barcode <- sub("-[0-9]+$", "", clust_df$barcode)
    
    # 4. CRITICAL FIX: Sync the spatial coordinate rownames with new barcodes
    rownames(image_coords) <- seu_clean
  } else {
    stop("Critical Error: Barcodes do not match even after cleaning.")
  }
}

# --- CREATE IMAGE OBJECT (Must be after barcode renaming) ---
colnames(image_coords) <- c("col", "row")
seu[['image']] <- new(Class = 'SlideSeq', 
                      assay = "Spatial", 
                      coordinates = as.data.frame(image_coords))

# 6. Run Annotation Loop
truth_cols <- setdiff(names(truth_df), c("X", "ID", "barcode"))

for (tc in truth_cols) {
  message(paste("Processing truth column:", tc))
  
  label_map <- truth_df[[tc]]
  names(label_map) <- truth_df$barcode
  groups <- label_map[colnames(seu)]
  
  valid_cells <- !is.na(groups)
  if (sum(valid_cells) == 0) next
  
  sub_expr <- expr_mat[, valid_cells, drop=FALSE]
  sub_groups <- groups[valid_cells]
  group_indices <- split(seq_len(ncol(sub_expr)), sub_groups)
  ref_avg <- do.call(cbind, lapply(group_indices, function(idx) rowMeans(sub_expr[, idx, drop=FALSE])))
  
  se_ref <- SummarizedExperiment(
    assays  = SimpleList(logcounts = as.matrix(ref_avg)),
    colData = DataFrame(label = colnames(ref_avg), row.names = colnames(ref_avg))
  )
  
  sr <- SingleR(test = expr_mat, ref = se_ref, labels = se_ref$label,
                assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  pred_df <- data.frame(barcode = colnames(expr_mat), sr_label = sr$labels, stringsAsFactors = FALSE)
  
  merged <- clust_df %>%
    inner_join(truth_df %>% select(barcode, !!sym(tc)), by = "barcode") %>%
    inner_join(pred_df, by = "barcode")
  
  if (nrow(merged) == 0) next
  
  annot <- merged %>%
    group_by(cluster = !!sym(tail(names(clust_df), 1))) %>%
    count(sr_label) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup()
  
  meta_col <- paste0("annotation_", tc)
  seu@meta.data[[meta_col]] <- factor(
    annot$sr_label[match(seu@meta.data$cluster, annot$cluster)]
  )
  
  # PREVENT CRASH: Ensure NAs are a named level so SpatialDimPlot/LabelClusters doesn't break
  seu@meta.data[[meta_col]] <- addNA(seu@meta.data[[meta_col]])
  levels(seu@meta.data[[meta_col]])[is.na(levels(seu@meta.data[[meta_col]]))] <- "Unassigned"
}

# 7. Save and Plot
if (tolower(opt$save_rds) %in% c("true", "yes", "1")) {
  dir.create(dirname(opt$seurat_out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seu, file = opt$seurat_out)
  message("Saved annotated Seurat object to: ", opt$seurat_out)
} else {
  dir.create(dirname(opt$seurat_out), showWarnings = FALSE, recursive = TRUE)
  writeLines("Annotation RDS saving skipped by config.", opt$seurat_out)
  message("Wrote placeholder to: ", opt$seurat_out, " (save_rds disabled)")
}

for (tc in truth_cols) {
  meta_col <- paste0('annotation_', tc)
  if (meta_col %in% colnames(seu@meta.data)) {
    # Check if there are valid categories to plot
    if (length(unique(na.omit(seu@meta.data[[meta_col]]))) > 0) {
      p <- SpatialDimPlot(seu, group.by = meta_col, label = TRUE) + 
        labs(title = paste('Cell Type Annotation (Ref:', tc, ')'))
      
      plot_path <- sub('\\.png$', paste0('_', tc, '.png'), opt$plot_seurat)
      dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)
      ggsave(plot_path, p, width = 16, height = 14)
    }
  }
}

message("Annotation step completed successfully.")