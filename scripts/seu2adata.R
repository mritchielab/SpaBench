library(reticulate)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(rhdf5)

assay <- DefaultAssay(seu)
VariableFeatures(seu, assay = assay) <- SpatiallyVariableFeatures(seu, assay= assay, selection.method=feat)

SaveH5Seurat(seu, filename = "/vast/projects/SpatialTranscriptomicsBenchmarking/Snakemake_workflow/data/HCC/test.h5Seurat")
Convert("/vast/projects/SpatialTranscriptomicsBenchmarking/Snakemake_workflow/data/HCC/test.h5Seurat", dest = "h5ad")

SaveH5Seurat(
  seu,
  filename     = "/vast/projects/SpatialTranscriptomicsBenchmarking/data/HCC/seu.h5Seurat",
  overwrite     = TRUE,
  write.data   = TRUE,
  assay        = c("Spatial"),
  write.images = TRUE
)

# 再转换
Convert(
  "/vast/projects/SpatialTranscriptomicsBenchmarking/test_sct.h5Seurat",
  assay = assay,
  dest      = "h5ad",
  overwrite = TRUE
)




# 1. 将 v5 Spatial assay 转成“旧版” Assay 对象
#    这里先在 seu 对象里新建一个 “Spatial3” assay
seu[["Spatial3"]] <- as(object = seu[["Spatial"]], Class = "Assay")

# 2. 把默认 assay 切换到 Spatial3，然后删掉旧的 Spatial
DefaultAssay(seu) <- "Spatial3"
seu[["Spatial"]] <- NULL

# 3. 把 “Spatial3” 重命名回 “Spatial” —— 这样内部是 v3 assay，但名字保持一致
seu <- RenameAssays(object = seu, Spatial3 = "Spatial")

# 检查一下
seu
# 应该看到 Active assay: Spatial (Assay v3)

# 4. 再用 SeuratDisk 把表达矩阵写入 .h5Seurat
SaveH5Seurat(
  seu,
  filename     = "/vast/projects/SpatialTranscriptomicsBenchmarking/data/HCC/test_v3Spatial.h5Seurat",
  overwrite     = TRUE,
  assay        = "Spatial",    # 现在这个 Spatial 已经是 v3 assay 了
  write.data   = TRUE,
  write.counts = TRUE,
  write.scale  = TRUE,
  write.images = TRUE
)

# 5. 最后用 Convert 导出 h5ad
Convert(
  source    = "/vast/projects/SpatialTranscriptomicsBenchmarking/test_v3Spatial.h5Seurat",
  dest      = "/vast/projects/SpatialTranscriptomicsBenchmarking/test_v3Spatial.h5ad",
  assay     = "Spatial",
  overwrite = TRUE,
  verbose   = TRUE
)






####
library(Seurat)
library(SeuratDisk)

# （1）读入原始 Seurat 对象
seu <- readRDS("/vast/projects/SpatialTranscriptomicsBenchmarking/Snakemake_workflow/data/Human_DLPFC_151510/01_raw_seurat.rds")

# （2）如果需要，先做归一化
seu <- NormalizeData(seu, assay="Spatial", verbose=FALSE)

# （3）手动降级为老版 Assay
Assay3 <- CreateAssayObject(counts = GetAssayData(seu, assay="Spatial", slot="counts"))
Assay3@data <- GetAssayData(seu, assay="Spatial", slot="data")

feature_df <- data.frame(
  gene_ids      = rownames(Assay3),
  feature_types = rep("Gene Expression", nrow(Assay3)),
  genome        = rep("hg38", nrow(Assay3)),
  row.names     = rownames(Assay3),
  stringsAsFactors = FALSE
)
Assay3@meta.features <- feature_df

seu[["Spatial3"]] <- Assay3
DefaultAssay(seu) <- "Spatial3"
seu[["Spatial"]] <- NULL
seu <- RenameAssays(seu, Spatial3 = "Spatial")

md <- seu@meta.data
md$array_row   <- coords[rownames(md), "row"]
md$array_col   <- coords[rownames(md), "col"]
md$in_tissue   <- TRUE
# 只保留这三列，并且顺序和 bdata 一致
seu@meta.data  <- md[, c("in_tissue", "array_row", "array_col")]

# （5）保存为 v3‐compatible h5Seurat，extra 把 images 和我们填好的 obs/var 一并存入
SaveH5Seurat(
  object       = seu,
  filename     = "/vast/projects/SpatialTranscriptomicsBenchmarking/data/HCC/out_v3Spatial.h5Seurat",
  overwrite    = TRUE,
  assay        = "Spatial",
  write.data   = TRUE,
  write.counts = TRUE,
  write.images = TRUE
)

# （6）再转成 h5ad
Convert(
  source    = "/vast/projects/SpatialTranscriptomicsBenchmarking/out_v3Spatial.h5Seurat",
  dest      = "/vast/projects/SpatialTranscriptomicsBenchmarking/out_v3Spatial.h5ad",
  assay     = "Spatial",
  overwrite = TRUE,
  verbose   = TRUE
)



#####


# （4）**新增：把 @images[[1]]@coordinates 写入 meta.data**  
coords <- seu@images[["slice1"]]@coordinates
# coords 应该包含 row, col 两列
seu@meta.data$array_row <- coords[rownames(seu@meta.data), "row"]
seu@meta.data$array_col <- coords[rownames(seu@meta.data), "col"]
seu@meta.data$orig.ident <- NULL
seu$nCount_Spatial3 <- NULL
seu$nFeature_Spatial3 <- NULL




Convert("/vast/projects/SpatialTranscriptomicsBenchmarking/Data/Mouse_embryo_stereoseq/E10.5_E1S1.MOSTA.h5ad", 
        dest = "h5seurat", overwrite = TRUE, verbose = TRUE, assay = "Spatial")
seu <- LoadH5Seurat("/vast/projects/SpatialTranscriptomicsBenchmarking/Data/Mouse_embryo_stereoseq/E10.5_E1S1.MOSTA.h5seurat", assays = "Spatial", image = FALSE)
barcodes <- rownames(seu@meta.data)
nums_list <- regmatches(barcodes, gregexpr("\\d+", barcodes))
coords <- as.data.frame(do.call(rbind, lapply(nums_list, function(v) {
  v <- as.numeric(v)
  if(length(v) >= 2) v[1:2] else c(NA, NA)
})))
colnames(coords) <- c("imagecol","imagerow")
coords$row <- coords$imagerow
coords$col <- coords$imagecol
rownames(coords) <- rownames(seu@meta.data)

seu[['slice1']] <- new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = coords)
SpatialDimPlot(seu)

