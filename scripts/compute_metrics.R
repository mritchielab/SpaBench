#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(aricode)
  library(mclust)
  library(FNN)
  library(cluster)
})

option_list <- list(
  make_option("--seurat_annot", type="character", help="Annotated Seurat RDS"),
  make_option("--truth",        type="character", help="Ground-truth TSV"),
  make_option("--out",          type="character", help="Metrics CSV output"),
  make_option("--plot",         type="character", help="Metrics plot PNG"),
  make_option("--algorithm",    type="character", help="Algorithm name")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load Data
seu <- readRDS(opt$seurat_annot)
truth_df <- read.delim(opt$truth, stringsAsFactors=FALSE)

# 2. Barcode Harmonization Function
clean_barcodes <- function(x) {
  x <- sub(" .*$", "", x)           # Remove " Cavity" etc.
  x <- sub("-[0-9]+$", "", x) # Remove the 10x "-1" suffix only (keeps Stereo-seq "row_col")
  return(x)
}

# Clean truth barcodes once
if ("ID" %in% names(truth_df)) {
  truth_df$barcode_clean <- clean_barcodes(as.character(truth_df$ID))
} else {
  truth_df$barcode_clean <- clean_barcodes(as.character(truth_df$X))
}

# 3. Define Metric Helper Functions (CHAOS, PAS, etc. - kept for brevity)
homogeneity <- function(true, pred) { tbl <- table(pred, true); sum(apply(tbl,1,max)) / length(true) }
completeness <- function(true, pred) { tbl <- table(pred, true); sum(apply(tbl,2,max)) / length(true) }
compute_chaos <- function(coords, labels, k = 1, normalize_max = NULL) {
  if (!is.null(normalize_max)) {
    coords <- apply(coords, 2, function(x) (x - min(x)) / diff(range(x)))
    coords <- coords * normalize_max
  }
  chaos <- 0
  for (cl in unique(labels)) {
    idx <- which(labels == cl); if (length(idx) > 1) {
      sub_coords <- coords[idx, , drop = FALSE]
      nn <- FNN::get.knn(sub_coords, k = k)
      chaos <- chaos + sum(nn$nn.dist)
    }
  }
  return(chaos)
}
compute_pas <- function(coords, labels, k = 10, threshold = 6) {
  nn_idx <- get.knn(coords, k = k)$nn.index
  diffs <- sapply(seq_len(nrow(nn_idx)), function(i) { sum(labels[nn_idx[i, ]] != labels[i]) })
  mean(diffs >= threshold)
}
compute_asw_safe <- function(coords, labels, max_total = 10000) {
  n <- nrow(coords)
  if (n > max_total) {
    set.seed(42)
    idx <- sample(seq_len(n), max_total)
    coords <- coords[idx,]; labels <- labels[idx]
  }
  dmat <- dist(coords)
  sil <- silhouette(as.numeric(factor(labels)), dmat)
  return((mean(sil[, "sil_width"], na.rm = TRUE) + 1) / 2)
}

# 4. Main Metric Loop
truth_cols <- setdiff(names(truth_df), c('X','ID','barcode','barcode_clean'))
metrics_list <- list()

for(tc in truth_cols) {
  annot_col <- paste0('annotation_', tc)
  
  if (!annot_col %in% colnames(seu@meta.data)) {
    message(paste("Skipping", tc, "because", annot_col, "is missing."))
    next
  }
  
  # --- ALIGNMENT WITHOUT RENAMING SEURAT ---
  # Create a temporary dataframe with original and clean barcodes
  df_map <- data.frame(
    barcode_orig = colnames(seu),
    barcode_clean = clean_barcodes(colnames(seu)),
    pred = seu@meta.data[[annot_col]],
    stringsAsFactors = FALSE
  )
  
  # Join with truth based on the cleaned barcodes
  merged <- df_map %>% 
    inner_join(truth_df, by = 'barcode_clean')
  
  if (nrow(merged) == 0) {
    message(paste("No overlapping barcodes for truth column:", tc))
    next
  }
  
  true <- merged[[tc]]
  pred <- merged$pred
  
  # --- COORDINATE EXTRACTION ---
  # annotate_clusters stores its coordinates as image "image"; use it when present
  img_name <- if ("image" %in% names(seu@images)) "image" else names(seu@images)[1]
  img_coords_all <- seu@images[[img_name]]@coordinates
  # Use the original barcode to pull coordinates, then align with merged data
  img_coords <- img_coords_all[merged$barcode_orig, ]
  
  available <- colnames(img_coords)
  if (all(c("imagerow", "imagecol") %in% available)) {
    selected_cols <- c("imagerow", "imagecol")
  } else if (all(c("row", "col") %in% available)) {
    selected_cols <- c("row", "col")
  } else {
    selected_cols <- colnames(img_coords)[1:2]
  }
  
  coords <- as.matrix(img_coords[, selected_cols])
  
  # Calculate Metrics
  metrics_list[[tc]] <- data.frame(
    truth=tc, 
    NMI=NMI(true,pred), 
    ARI=adjustedRandIndex(true,pred),
    HOM=homogeneity(true,pred), 
    COM=completeness(true,pred),
    CHAOS=compute_chaos(coords, pred, normalize_max=100),
    PAS = compute_pas(coords, pred),
    ASW = tryCatch(compute_asw_safe(coords, pred), error = function(e) NA)
  )
}

# 5. Output
metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, opt$out, row.names=FALSE)

# Plotting
library(ggplot2)
long <- tidyr::pivot_longer(metrics_df, cols=c('NMI','ARI','HOM','COM'), names_to='Metric', values_to='Value')
p <- ggplot(long, aes(x=truth, y=Value, fill=Metric)) + 
  geom_col(position='dodge') + facet_wrap(~Metric, scales='free_y') + 
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title=paste('Metrics:', opt$algorithm))

ggsave(opt$plot, p, width=12, height=10)