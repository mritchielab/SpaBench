#!/usr/bin/env Rscript
#
# aggregate_cluster_metrics.R
#
# Concatenates every results/<sample>/<norm>/<feat>/<alg>/09_<alg>_cluster_metrics.csv
# into a single results/<sample>/clustering_metrics.csv with explicit
# Sample / Normalization / Feature_selection / Clustering identifier columns.
#
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--metrics_dir", type = "character",
              help = "Sample directory, e.g. results/Mouse_brain_squidpy"),
  make_option("--out_csv",     type = "character", help = "Output CSV path")
)))

if (is.null(opt$metrics_dir) || is.null(opt$out_csv)) {
  stop("--metrics_dir and --out_csv are required")
}

files <- list.files(
  path       = opt$metrics_dir,
  pattern    = "^09_.*_cluster_metrics\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("No 09_*_cluster_metrics.csv files found under ", opt$metrics_dir)

# Sample is the basename of metrics_dir.
sample_name <- basename(normalizePath(opt$metrics_dir))

df_list <- lapply(files, function(f) {
  # Path layout: <metrics_dir>/<norm>/<feat>/<alg>/09_<alg>_cluster_metrics.csv
  rel <- sub(paste0("^", normalizePath(opt$metrics_dir), "/"), "",
             normalizePath(f))
  parts <- strsplit(rel, "/")[[1]]
  if (length(parts) < 4) return(NULL)
  norm_v  <- parts[1]
  feat_v  <- parts[2]
  alg_v   <- parts[3]
  body <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(body) == 0) return(NULL)
  data.frame(
    Sample            = sample_name,
    Normalization     = norm_v,
    Feature_selection = feat_v,
    Clustering        = alg_v,
    body,
    stringsAsFactors  = FALSE
  )
})
all <- bind_rows(df_list)

dir.create(dirname(opt$out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(all, opt$out_csv, row.names = FALSE)

cat(sprintf("Wrote %s  (rows = %d, pipelines = %d)\n",
            opt$out_csv, nrow(all),
            length(unique(paste(all$Normalization, all$Feature_selection, all$Clustering)))))
