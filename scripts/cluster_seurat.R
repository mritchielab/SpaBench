# scripts/cluster_seurat.R
#!/usr/bin/env Rscript
orig_options <- base::options; options <- function(...){dots<-list(...); if(!is.null(dots$scipen)){max_sci<-.Machine$integer.max-1; dots$scipen<-min(as.integer(dots$scipen),max_sci)}; do.call(orig_options,dots)}
Sys.setenv(RETICULATE_PYTHON="/stornext/System/data/software/rhel/9/base/tools/python/3.11.8/bin/python")
suppressMessages(library(optparse)); library(Seurat); library(patchwork); library(ggplot2)
opt<-parse_args(OptionParser(option_list=list(
  make_option("--input",      type="character"),
  make_option("--output",     type="character"),
  make_option("--resolution", type="numeric"),
  make_option("--algorithm",  type="integer"),
  make_option("--clustermethod",type="character",default="matrix"),
  make_option("--meta_csv",   type="character"),
  make_option("--plot_png",   type="character")
)))
seu<-readRDS(opt$input)

# features used for clustering
feat <- basename(dirname(opt$input))
message("Feature selection = ", feat)
assay <- DefaultAssay(seu)

if (!(feat %in% c("vst", "nnSVG", "SPARK", "SPARK-X", "disp", "mvp"))) {
  # moransi / markvariogram
  if (feat == "moransi") {
    message("Overriding var.features with SpatiallyVariableFeatures(method = 'moransi')")
    VariableFeatures(seu, assay = assay) <- SpatiallyVariableFeatures(
      seu, assay = assay, selection.method = "moransi"
    )
  } else if (feat == "markvariogram") {
    message("Overriding var.features with SpatiallyVariableFeatures(method = 'markvariogram')")
    VariableFeatures(seu, assay = assay) <- SpatiallyVariableFeatures(
      seu, assay = assay, selection.method = "markvariogram"
    )
  } else {
    stop("Unsupported feature‐selection method: ", feat)
  }
  
} else if (feat %in% c("vst", "disp", "mvp")) {
  message("Using existing VariableFeatures (vst / disp / mvp)")
  
} else if (feat %in% c("nnSVG", "SPARK", "SPARK-X")) {
  message(sprintf(
    "Using features from %s (already set in VariableFeatures(seu))",
    feat
  ))
  
} else {
  stop("Unknown feature‐selection method: ", feat)
}

print(opt$algorithm)
npcs<-ncol(Embeddings(seu,"pca"))
seu<-FindNeighbors(seu,dims=1:npcs)
seu<-FindClusters(seu,resolution=opt$resolution,algorithm=opt$algorithm,method=opt$clustermethod)
write.csv(seu@meta.data,opt$meta_csv,row.names=TRUE)
saveRDS(seu,opt$output)
p1<-DimPlot(seu,reduction="umap",label=TRUE)
p2<-SpatialDimPlot(seu,label=TRUE,label.size=4.5)
ggsave(opt$plot_png,p1+p2,width=14,height=12)
