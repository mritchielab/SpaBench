# scripts/aggregate_metrics.R

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--metrics_dir", type="character", help="根目录（results/{sample}）"),
  make_option("--out_csv",      type="character", help="Output CSV path"),
  make_option("--out_plot",     type="character", help="Output PNG oath")
)))

# 1. 找到所有单算法 metrics 文件
files <- list.files(
  path       = file.path(opt$metrics_dir),
  pattern    = "^09_.*_metrics\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)
# The pattern also matches 09_<alg>_cluster_metrics.csv (main route); keep the SingleR files only
files <- files[!grepl("_cluster_metrics\\.csv$", files)]

# 2. 读取并打上标签（从路径中提取 norm/feat/red/alg）
df_list <- lapply(files, function(f) {
  parts <- strsplit(f, "/")[[1]]
  # parts: ".../results/{sample}/{norm}/{feat}/{alg}/09_{alg}_metrics.csv"
  n <- tail(parts, 4)
  tibble(
    norm   = n[1],
    feat   = n[2],
    alg    = n[3],
    read.csv(f, stringsAsFactors=FALSE)
  )
})
all_metrics <- bind_rows(df_list)
write.csv(all_metrics, opt$out_csv, row.names=FALSE)

long_df <- all_metrics %>%
  pivot_longer(cols = c("NMI","ARI","HOM","COM","CHAOS","PAS","ASW"),
               names_to  = "Metric",
               values_to = "Value") %>%
  unite("Pipeline", norm:feat:red, sep="|", remove=FALSE)

p <- ggplot(long_df, aes(x=Pipeline, y=Value, fill=alg)) +
  geom_col(position="dodge") +
  facet_wrap(~Metric, scales="free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "bottom"
  ) +
  labs(x="Pipeline (norm|feat|red)", y="Score", title="All clustering metrics")

ggsave(opt$out_plot, p, width=16, height=14)

