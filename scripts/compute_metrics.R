#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(aricode)  # NMI
  library(mclust)   # ARI
  library(FNN)      # get.knn
  library(cluster)  # silhouette
})

option_list <- list(
  make_option("--seurat",      type="character", help="Input Seurat RDS (04_reduced.rds)"),
  make_option("--clusters",    type="character", help="Clustering CSV (07_*_clusters.csv)"),
  make_option("--truth",       type="character", help="Ground‑truth TSV"),
  make_option("--out",         type="character", help="Metrics CSV output (09_*_metrics.csv)"),
  make_option("--plot",        type="character", default=NULL, help="(Optional) single‐task metrics plot PNG"),
  make_option("--seurat_out",  type="character", help="Output updated Seurat RDS with annotation"),
  make_option("--plot_seurat", type="character", help="Output PNG for Seurat annotation plot")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. read in Seurat
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
  image_coords <- seu@meta.data[, c("imagecol","imagerow"), drop = FALSE]
  colnames(image_coords) <- c("col", "row")
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
}


# 3. Load ground truth
truth_df <- read.delim(opt$truth, stringsAsFactors=FALSE)
if ("ID" %in% names(truth_df)) {
  truth_df$barcode <- truth_df$ID
} else if ("X" %in% names(truth_df)) {
  truth_df$barcode <- truth_df$X
} else {
  stop("Cannot find ID or X column in truth file")
}

# 4. Metrics helpers
homogeneity <- function(true, pred) {
  tbl <- table(pred, true)
  sum(apply(tbl,1,max)) / length(true)
}
completeness <- function(true, pred) {
  tbl <- table(pred, true)
  sum(apply(tbl,2,max)) / length(true)
}
# CHAOS：对每个 cluster，计算同一 cluster 内所有 spot 两两距离之和，最后加总
compute_chaos <- function(coords, labels, k = 1) {
  N <- length(labels)
  total_dist <- 0
  
  for (cl in unique(labels)) {
    idx <- which(labels == cl)
    if (length(idx) <= k) next
    
    sub_coords <- coords[idx, , drop = FALSE]
    # get.knn 返回最近 k 个邻居的距离矩阵
    nn <- get.knn(sub_coords, k = k)
    # nn$nn.dist 是一个 m × k 矩阵，m = nrow(sub_coords)
    # 把所有这些距离加起来
    total_dist <- total_dist + sum(nn$nn.dist)
  }
  
  return(total_dist / N)
}

# PAS：对每个 spot 看它 10 最近邻中，与它标签不同的数量是否 ≥6
compute_pas <- function(coords, labels, k = 10, threshold = 6) {
  # 得到一个 n × k 的矩阵 nn_idx，行号是 spot，内容是它的 k 个邻居
  nn_idx <- get.knn(coords, k = k)$nn.index
  
  # 对每一行 i，计算 labels[nn_idx[i, ]] != labels[i] 的个数
  diffs <- sapply(seq_len(nrow(nn_idx)), function(i) {
    sum(labels[nn_idx[i, ]] != labels[i])
  })
  
  # 最后返回“不相等个数 ≥ threshold 的 spot 百分比”
  mean(diffs >= threshold)
}

# ASW：用物理距离矩阵和预测标签算 silhouette，归一化到 [0,1]
compute_asw <- function(coords, labels) {
  dmat <- dist(coords)
  sil <- silhouette(as.numeric(factor(labels)), dmat)
  # silhouette 带负值，这里我们取平均并重新缩放到 0–1
  mean(sil[, "sil_width"] + 1) / 2
}

# 5. Compute metrics per truth column
truth_cols   <- setdiff(names(truth_df), c("X","ID","barcode"))
cluster_maps <- list()
debug_sr     <- NULL   
debug_merged <- NULL   

metrics_list <- lapply(truth_cols, function(tc) {
  # 5.1 build pseudo‑bulk reference
  groups  <- truth_df[[tc]]
  ref_avg <- do.call(cbind, lapply(split(seq_len(ncol(expr_mat)), groups), function(idx) {
    rowMeans(expr_mat[, idx, drop=FALSE])
  }))
  colnames(ref_avg) <- names(split(seq_len(ncol(expr_mat)), groups))
  rownames(ref_avg) <- rownames(expr_mat)
  se_ref <- SummarizedExperiment(
    assays  = SimpleList(logcounts=ref_avg),
    colData = DataFrame(label=colnames(ref_avg), row.names=colnames(ref_avg))
  )
  
  # 5.2 run SingleR
  sr <- SingleR(
    test            = expr_mat,
    ref             = se_ref,
    labels          = se_ref$label,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts"
  )
  debug_sr <<- sr             ## <<< MOD
  warning(sprintf("Finished SingleR for truth column '%s', top labels: %s", 
                  tc, paste(head(sr$labels), collapse=",")))
  
  pred_df <- data.frame(barcode=colnames(expr_mat),
                        sr_label=sr$labels,
                        stringsAsFactors=FALSE)
  
  # 5.3 merge
  merged <- clust_df %>%
    inner_join(truth_df, by="barcode") %>%
    inner_join(pred_df,  by="barcode")
  debug_merged <<- merged     ## <<< MOD
  
  # 5.4 cluster annotation map
  cluster_col <- tail(names(clust_df),1)
  annot <- merged %>%
    group_by_at(cluster_col) %>%
    count(sr_label) %>%
    slice_max(n, with_ties=FALSE) %>%
    select(!!sym(cluster_col), annot=sr_label)
  # <<< MOD: 用 <<- 保证给外部的列表赋值
  cluster_maps[[tc]] <<- annot 
  
  merged <- merged %>% left_join(annot, by=cluster_col)
  
  # 5.5 compute metrics
  #coords <- as.matrix(seu@images[[1]]@coordinates[merged$barcode, c("row","col")])
  if (any(c("pred", "refined_pred") %in% colnames(seu@meta.data))) {
        coords <- as.matrix(
            seu@meta.data[ merged$barcode, c("x_pixel", "y_pixel") ]
          )
       } else if ("stlearn" %in% colnames(seu@meta.data)) {
          coords <- as.matrix(
            # 
            seu@meta.data[ merged$barcode, c("array_row", "array_col") ]
          )
      } else {
          # normal workflows within Seurat then use images slot
            coords <- as.matrix(seu@images[[1]]@coordinates[ merged$barcode, c("row", "col") ])
          }
  pred   <- merged$annot
  chaos_val <- compute_chaos(coords, pred)
  pas_val   <- compute_pas(coords, pred)
  asw_val   <- compute_asw(coords, pred)
  
  true_vec <- merged[[tc]]
  pred_vec <- merged$annot
  keep     <- true_vec!="unknown"
  data.frame(
    truth = tc,
    NMI   = NMI(true_vec[keep], pred_vec[keep]),
    ARI   = adjustedRandIndex(true_vec[keep], pred_vec[keep]),
    HOM   = homogeneity(true_vec[keep], pred_vec[keep]),
    COM   = completeness(true_vec[keep], pred_vec[keep]),
    CHAOS = chaos_val,
    PAS   = pas_val,
    ASW   = asw_val,
    stringsAsFactors=FALSE
  )
})

metrics_df <- bind_rows(metrics_list)
write.csv(metrics_df, opt$out, row.names=FALSE)

# 6. Optional single‐task metrics plot
if (!is.null(opt$plot)) {
  long_df <- pivot_longer(metrics_df,
                          cols=c("NMI","ARI","HOM","COM"),
                          names_to="Metric", values_to="Value")
  p <- ggplot(long_df, aes(x=truth, y=Value, fill=Metric)) +
    geom_col(position="dodge") +
    facet_wrap(~Metric, scales="free_y") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="none") +
    labs(x="Ground truth column", y="Score",
         title="Clustering metrics by truth")
  ggsave(opt$plot, p, width=16, height=14)
}

# 7. Annotate Seurat and plot for each truth column
cluster_col <- tail(names(clust_df), 1)

dir.create(dirname(opt$seurat_out), showWarnings = FALSE, recursive = TRUE)
plots <- list()
for (tc in truth_cols) {
  # truth corresponding cluster->annotation map
  annot_map <- cluster_maps[[tc]]
  # add metadata as annotation_<tc>
  meta_col <- paste0("annotation_", tc)
  seu@meta.data[[meta_col]] <- factor(
    annot_map$annot[ match(
      as.character(seu@meta.data[[cluster_col]]),
      annot_map[[cluster_col]]
    )]
  )
  p <- SpatialDimPlot(seu, group.by = meta_col, label = TRUE) +
    labs(title = paste("Annotation by", tc))
  plot_file <- sub("\\.rds$", paste0("_", tc, ".png"), opt$plot_seurat)
  ggsave(plot_file, p, width = 16, height = 14)
  plots[[tc]] <- p
}
saveRDS(seu, file=opt$seurat_out)

# 8. write out debug objects for interactive inspection
saveRDS(debug_sr,     file="debug_last_sr.rds")     ## <<< MOD
saveRDS(debug_merged, file="debug_last_merged.rds") ## <<< MOD
