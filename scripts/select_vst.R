# scripts/select_vst.R
#!/usr/bin/env Rscript
# Highly variable feature selection with Seurat's vst method (top n_features by
# standardized variance).
#  * Standard assays: FindVariableFeatures(selection.method = "vst") (counts layer).
#  * SCTAssay: FindVariableFeatures() would dispatch to the SCT method and ignore
#    selection.method, so vst is computed explicitly on the SCT corrected counts
#    (SCT counts layer) and set as VariableFeatures of the SCT assay.
#  * No layer of the input object is modified.
suppressMessages(library(optparse)); library(Seurat); library(ggplot2)
option_list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--method", type="character"),
  make_option("--n_features", type="integer"),
  make_option("--csv", type="character"),
  make_option("--plot", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
seu <- readRDS(opt$input)
assay <- DefaultAssay(seu)
is_sct <- inherits(seu[[assay]], "SCTAssay")
message("Assay: ", assay, " (class ", class(seu[[assay]])[1], ")")

if (is_sct) {
  message("SCTAssay detected: computing vst explicitly on SCT corrected counts")
  counts <- LayerData(seu, assay = assay, layer = "counts")
  if (!inherits(counts, "dgCMatrix")) counts <- as(as.matrix(counts), "CsparseMatrix")
  hvf <- Seurat:::VST(data = counts, nselect = opt$n_features)
  rownames(hvf) <- rownames(counts)
  rm(counts)
  ok <- !is.na(hvf$rank)
  selected <- rownames(hvf)[ok][order(hvf$rank[ok])]
  VariableFeatures(seu[[assay]]) <- selected
} else {
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=opt$n_features)
  selected <- VariableFeatures(seu)
}
stopifnot(identical(VariableFeatures(seu), selected))
message("Selected ", length(selected), " features (requested ", opt$n_features, ")")

write.csv(selected, opt$csv, row.names=FALSE)
saveRDS(seu, opt$output)

# plot variable features
if (is_sct) {
  hvf$selected <- rownames(hvf) %in% selected
  p <- ggplot(hvf, aes(x = mean, y = variance.standardized, colour = selected)) +
    geom_point(size = 0.5) + scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "red")) +
    labs(x = "Average corrected count", y = "Standardized variance",
         title = sprintf("vst on SCT corrected counts: %d selected features", length(selected))) +
    theme_bw()
} else {
  p <- VariableFeaturePlot(seu)
}
ggsave(opt$plot, p)
