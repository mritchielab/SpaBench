#!/usr/bin/env Rscript
#
# compute_cluster_metrics.R
#
# Computes the same accuracy / spatial-coherence metrics as compute_metrics.R,
# but uses the *raw cluster IDs* from 07_<alg>_clusters.csv as the prediction
# (instead of the SingleR-derived annotation_<truth> column, which leaks the
# ground-truth labels into the prediction).
#
# Inputs
#   --clusters    07_<alg>_clusters.csv  (first col = barcode, last col = cluster)
#   --positions   data/<sample>/01_positions.csv  (no header; col1 = barcode,
#                                                  cols 5-6 = pixel x/y)
#   --truth       data/Ground_truth/<sample>/labels.tsv
#   --out         output metrics CSV
#   --plot        output plot PNG
#   --algorithm   algorithm name (used only in plot title)
#
# Output CSV schema (one row per truth column, identical to compute_metrics.R):
#   truth, NMI, ARI, HOM, COM, CHAOS, PAS, ASW
#
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(aricode)
  library(mclust)
  library(FNN)
  library(cluster)
})

option_list <- list(
  make_option("--clusters",  type = "character", help = "07_*_clusters.csv"),
  make_option("--positions", type = "character", help = "01_positions.csv"),
  make_option("--truth",     type = "character", help = "Ground-truth TSV/CSV"),
  make_option("--out",       type = "character", help = "Output metrics CSV"),
  make_option("--plot",      type = "character", help = "Output plot PNG"),
  make_option("--algorithm", type = "character", default = "unknown",
              help = "Algorithm name (for plot title only)"),
  make_option("--ground_truth_label", type = "character", default = "full",
              help = paste("'full' (default) – every clustered spot has a truth",
                           "label and only labelled spots feed every metric.",
                           "'partial' – only some spots are labelled; missing",
                           "or empty labels are filled with '_Other_' so",
                           "NMI/ARI/HOM/COM stay well defined on the full",
                           "clustered set; CHAOS/PAS/ASW use the full set",
                           "in both modes."))
)
opt <- parse_args(OptionParser(option_list = option_list))

# ------- helpers ----------------------------------------------------------
clean_barcodes <- function(x) {
  x <- sub(" .*$", "", x)
  # Only strip the 10x GEM suffix ("-1"). Stripping "_<digit>" merged distinct
  # Stereo-seq barcodes such as "373_8" and "373_9".
  x <- sub("-[0-9]+$", "", x)
  x
}

homogeneity  <- function(true, pred) { tbl <- table(pred, true); sum(apply(tbl, 1, max)) / length(true) }
completeness <- function(true, pred) { tbl <- table(pred, true); sum(apply(tbl, 2, max)) / length(true) }

compute_chaos <- function(coords, labels, k = 1, normalize_max = NULL) {
  if (!is.null(normalize_max)) {
    coords <- apply(coords, 2, function(x) (x - min(x)) / diff(range(x)))
    coords <- coords * normalize_max
  }
  chaos <- 0
  for (cl in unique(labels)) {
    idx <- which(labels == cl)
    if (length(idx) > 1) {
      sub_coords <- coords[idx, , drop = FALSE]
      nn <- FNN::get.knn(sub_coords, k = k)
      chaos <- chaos + sum(nn$nn.dist)
    }
  }
  chaos
}

compute_pas <- function(coords, labels, k = 10, threshold = 6) {
  nn_idx <- get.knn(coords, k = k)$nn.index
  diffs <- vapply(seq_len(nrow(nn_idx)),
                  function(i) sum(labels[nn_idx[i, ]] != labels[i]),
                  integer(1))
  mean(diffs >= threshold)
}

compute_asw_safe <- function(coords, labels, max_total = 10000) {
  n <- nrow(coords)
  if (n > max_total) {
    set.seed(42)
    idx <- sample(seq_len(n), max_total)
    coords <- coords[idx, ]; labels <- labels[idx]
  }
  if (length(unique(labels)) < 2) return(NA_real_)
  dmat <- dist(coords)
  sil <- silhouette(as.numeric(factor(labels)), dmat)
  (mean(sil[, "sil_width"], na.rm = TRUE) + 1) / 2
}

# ------- 1. Load clusters --------------------------------------------------
clu <- read.csv(opt$clusters, stringsAsFactors = FALSE,
                check.names = FALSE)
# First column is barcode (may have an empty header in some files); last is cluster
names(clu)[1] <- "barcode"
clu$cluster <- as.character(clu[[ncol(clu)]])
clu$barcode <- as.character(clu$barcode)
clu$bc      <- clean_barcodes(clu$barcode)

# ------- 2. Load truth -----------------------------------------------------
# Auto-detect tab vs comma delimiter from the header line so we accept either
# .tsv (e.g. labels.tsv) or .csv (e.g. Eye_only.csv) ground-truth files.
.truth_first <- readLines(opt$truth, n = 1L, warn = FALSE)
.truth_sep   <- if (grepl("\t", .truth_first)) "\t" else ","
truth <- read.table(opt$truth, sep = .truth_sep, header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE,
                    quote = "\"", fill = TRUE, na.strings = c("", "NA"))
# A header that starts with an empty field (e.g. DLPFC "\tlabel") gives an
# unnamed first column; name it "X" so it is recognised as the barcode column.
names(truth)[!nzchar(names(truth))] <- "X"
bc_col <- if ("ID" %in% names(truth)) "ID" else if ("X" %in% names(truth)) "X" else names(truth)[1]
truth$bc <- clean_barcodes(as.character(truth[[bc_col]]))

# ------- 3. Load positions -------------------------------------------------
# Preferred: 01_positions.csv (no header, col1 = barcode, last 2 cols = x/y).
# Fallback (Stereo-seq): barcodes themselves encode "<row>_<col>" coordinates.
pos <- NULL
if (!is.null(opt$positions) && nzchar(opt$positions) && file.exists(opt$positions)) {
  pos <- tryCatch(read.csv(opt$positions, stringsAsFactors = FALSE,
                           header = FALSE),
                  error = function(e) NULL)
}
if (!is.null(pos) && nrow(pos) > 0) {
  names(pos)[1] <- "barcode"
  pos$bc <- clean_barcodes(as.character(pos$barcode))
  nc <- ncol(pos) - 1
  pos$x <- suppressWarnings(as.numeric(pos[[nc - 1]]))
  pos$y <- suppressWarnings(as.numeric(pos[[nc]]))
  pos <- pos[, c("bc", "x", "y")]
} else {
  bc_sample <- head(clu$bc, 5)
  if (all(grepl("^[0-9]+_[0-9]+$", bc_sample))) {
    message("No positions.csv; parsing coordinates from Stereo-seq \"row_col\" barcode format.")
    pos <- data.frame(
      bc = clu$bc,
      x  = as.numeric(sub("_.*$", "", clu$bc)),
      y  = as.numeric(sub("^[^_]*_", "", clu$bc)),
      stringsAsFactors = FALSE
    )
  } else if (all(grepl("^s_[0-9]+um_[0-9]+_[0-9]+$", bc_sample))) {
    message("No positions.csv; parsing coordinates from VisiumHD \"s_NNNum_row_col\" barcode format.")
    parts <- strsplit(clu$bc, "_", fixed = TRUE)
    pos <- data.frame(
      bc = clu$bc,
      x  = as.numeric(vapply(parts, `[`, character(1), 3)),
      y  = as.numeric(vapply(parts, `[`, character(1), 4)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("No positions file and barcodes do not encode coordinates: ", opt$positions)
  }
}

# ------- 4. Build full clustered set (no truth filter yet) ----------------
full_df <- clu %>% select(bc, cluster) %>%
  inner_join(pos, by = "bc")
if (nrow(full_df) == 0) stop("No barcode overlap between clusters and positions.")

coords_full <- as.matrix(full_df[, c("x", "y")])
pred_full   <- full_df$cluster

# ------- 5. Compute metrics for every truth column ------------------------
truth_cols <- setdiff(names(truth), c("X", "ID", "barcode", "bc", bc_col))
metrics_list <- list()

partial_mode <- identical(tolower(opt$ground_truth_label), "partial")

for (tc in truth_cols) {
  truth_map <- setNames(as.character(truth[[tc]]), truth$bc)
  raw_truth <- truth_map[full_df$bc]   # NA for spots not in the truth file
  raw_truth[!is.na(raw_truth) & raw_truth == ""] <- NA  # blanks count as missing
  n_labeled <- sum(!is.na(raw_truth))
  n_total   <- nrow(full_df)

  if (partial_mode) {
    # Fill unknown/blank truth with "_Other_" so NMI/ARI are well-defined on
    # the full clustered set; CHAOS/PAS/ASW use the same full set.
    true_v   <- raw_truth
    true_v[is.na(true_v)] <- "_Other_"
    pred_v   <- pred_full
    coords_v <- coords_full
  } else {
    # Default 'full' mode — restrict every metric to labelled spots only,
    # matching the original behaviour.
    ok       <- !is.na(raw_truth)
    true_v   <- raw_truth[ok]
    pred_v   <- pred_full[ok]
    coords_v <- coords_full[ok, , drop = FALSE]
  }

  if (length(true_v) < 2 ||
      length(unique(true_v)) < 2 || length(unique(pred_v)) < 2) {
    message("Skipping truth=", tc, " (single class or empty after filter)")
    next
  }

  metrics_list[[tc]] <- data.frame(
    truth     = tc,
    NMI       = aricode::NMI(true_v, pred_v),
    ARI       = mclust::adjustedRandIndex(true_v, pred_v),
    HOM       = homogeneity(true_v, pred_v),
    COM       = completeness(true_v, pred_v),
    CHAOS     = compute_chaos(coords_v, pred_v, normalize_max = 100),
    PAS       = compute_pas(coords_v, pred_v),
    ASW       = tryCatch(compute_asw_safe(coords_v, pred_v),
                         error = function(e) NA_real_),
    n_labeled = n_labeled,
    n_total   = n_total
  )
}

if (length(metrics_list) == 0) stop("No metrics computed; aborting.")

metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, opt$out, row.names = FALSE)

# ------- 6. Plot -----------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
long <- tidyr::pivot_longer(metrics_df,
                            cols = c("NMI", "ARI", "HOM", "COM"),
                            names_to = "Metric", values_to = "Value")
p <- ggplot(long, aes(x = truth, y = Value, fill = Metric)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Cluster-only metrics:", opt$algorithm))
ggsave(opt$plot, p, width = 10, height = 8, dpi = 200)

cat(sprintf("Wrote %s (rows = %d, n_cells_used = %d)\n",
            opt$out, nrow(metrics_df), nrow(full_df)))
